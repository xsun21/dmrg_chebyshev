function mpsB = reduceD(mpsA, mpoX, DB, precision)
% Variationally compress an MPS/MPO to a smaller bond dimension.
%
%   INPUTS
%     mpsA      - Cell array {1,N}: the MPS (or MPO reshaped as MPS) to compress.
%     mpoX      - Cell array {1,N}: optional MPO to apply before compression.
%     DB        - Target bond dimension of the compressed MPS mpsB.
%     precision - Convergence criterion.
%
%   OUTPUT
%     mpsB - Cell array {1,N} of compressed MPS tensors with bond dim <= DB.
%

N = length(mpsA);
d = size(mpsA{1}, 3);   % physical dimension (may be d^2 for MPO-as-MPS)

% Initialize mpsB as a random right-canonical MPS of bond dimension DB
mpsB = createrandommps(N, DB, d);
mpsB = prepare(mpsB);   % right-canonical form required before left sweep

% Build the initial right-environment storage for the overlap <mpsB | mpoX | mpsA>
Cstorage = initCstorage(mpsB, mpoX, mpsA, N);

while 1
    Kvalues = [];

    % ------ Left sweep: sites 1 -> N-1 ------
    for j = 1:(N - 1)
        Cleft  = Cstorage{j};
        Cright = Cstorage{j+1};
        A = mpsA{j};
        if isempty(mpoX)
            X = reshape(eye(d), [1, 1, d, d]);   % identity MPO
        else
            X = mpoX{j};
        end

        % Solve single-site problem: optimal B{j} given fixed environments
        [B, K] = reduceD2_onesite(A, X, Cleft, Cright);

        % Left-canonicalize B and push gauge to right
        [B, U] = prepare_onesite(B, 'lr');
        mpsB{j} = B;
        Kvalues = [Kvalues, K];

        % Update left environment to include site j
        Cstorage{j+1} = updateCleft(Cleft, B, X, A);
    end

    % ------ Right sweep: sites N -> 2 ------
    for j = N:(-1):2
        Cleft  = Cstorage{j};
        Cright = Cstorage{j+1};
        A = mpsA{j};
        if isempty(mpoX)
            X = reshape(eye(d), [1, 1, d, d]);
        else
            X = mpoX{j};
        end

        % Solve single-site problem
        [B, K] = reduceD2_onesite(A, X, Cleft, Cright);

        % Right-canonicalize B and push gauge to left
        [B, U] = prepare_onesite(B, 'rl');
        mpsB{j} = B;
        Kvalues = [Kvalues, K];

        % Update right environment to include site j
        Cstorage{j} = updateCright(Cright, B, X, A);
    end

    % Convergence check 
    if std(Kvalues) / abs(mean(Kvalues)) < precision
        % Absorb final gauge into site 1 to complete normalization
        mpsB{1} = contracttensors(mpsB{1}, 3, 2, U, 2, 1);
        mpsB{1} = permute(mpsB{1}, [1, 3, 2]);
        break;
    end
end


% ==========================================================================
function [B, K] = reduceD2_onesite(A, X, Cleft, Cright)
% Solve the single-site compression problem analytically.

% Contract left environment with ket A (left bond of A contracted into Cleft)
Cleft = contracttensors(Cleft, 3, 3, A, 3, 1);
% Contract MPO X: physical ket index (4) and left MPO bond (1) contracted
Cleft = contracttensors(Cleft, 4, [2, 4], X, 4, [1, 4]);
% Contract right environment: MPO+ket bonds contracted with Cright
B = contracttensors(Cleft, 4, [3, 2], Cright, 3, [2, 3]);
B = permute(B, [1, 3, 2]);   % reorder to [Dl_B, Dr_B, d]

% K = -||B||^2 serves as a proxy for the compression quality
b = reshape(B, [prod(size(B)), 1]);
K = -b' * b;
