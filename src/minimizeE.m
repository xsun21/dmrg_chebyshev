function [E, mps] = minimizeE(hset, D, precision, mpsB)
% DMRG ground-state energy search (single-site, 1D).
%
%   INPUTS
%     hset      - (M x N) cell array of local operators. 
%     D         - Bond dimension.
%     precision - Convergence tolerance.
%     mpsB      - MPS to orthogonalize against (for excited states). 
%
%   OUTPUTS
%     E   - Ground state energy estimate.
%     mps - Cell array {1,N} of MPS tensors for the ground state.
%

[M, N] = size(hset);
d      = size(hset{1, 1}, 1);

% Initialize random MPS in right-canonical form
mps = createrandommps(N, D, d);
mps = prepare(mps);

% Build initial right-environment storage for <mps|H|mps>
Hstorage = initHstorage(mps, hset, d);
if ~isempty(mpsB)
    Cstorage = initCstorage(mps, [], mpsB, N);
end
P = [];   % orthogonal projector (empty = unconstrained)

while 1
    Evalues = [];

    % ------ Cycle 1: left sweep (sites 1 -> N-1) ------
    for j = 1:(N - 1)
        % Compute projector onto complement of mpsB at this site
        if ~isempty(mpsB)
            B      = mpsB{j};
            Cleft  = Cstorage{j};
            Cright = Cstorage{j+1};
            P      = calcprojector_onesite(B, Cleft, Cright);
        end

        Hleft  = Hstorage(:, j);
        Hright = Hstorage(:, j+1);
        hsetj  = hset(:, j);

        % Single-site energy minimization via sparse diagonalization
        [A, E] = minimizeE_onesite(hsetj, Hleft, Hright, P);

        % Left-canonicalize and push gauge right
        [A, U]  = prepare_onesite(A, 'lr');
        mps{j}  = A;
        Evalues = [Evalues, E];

        % Update left environments
        for m = 1:M
            h = reshape(hset{m, j}, [1, 1, d, d]);
            Hstorage{m, j+1} = updateCleft(Hleft{m}, A, h, A);
        end
        if ~isempty(mpsB)
            Cstorage{j+1} = updateCleft(Cleft, A, [], B);
        end
    end

    % ------ Cycle 2: right sweep (sites N -> 2) ------
    for j = N:(-1):2
        if ~isempty(mpsB)
            B      = mpsB{j};
            Cleft  = Cstorage{j};
            Cright = Cstorage{j+1};
            P      = calcprojector_onesite(B, Cleft, Cright);
        end

        Hleft  = Hstorage(:, j);
        Hright = Hstorage(:, j+1);
        hsetj  = hset(:, j);

        [A, E] = minimizeE_onesite(hsetj, Hleft, Hright, P);

        [A, U]  = prepare_onesite(A, 'rl');
        mps{j}  = A;
        Evalues = [Evalues, E];

        for m = 1:M
            h = reshape(hset{m, j}, [1, 1, d, d]);
            Hstorage{m, j} = updateCright(Hright{m}, A, h, A);
        end
        if ~isempty(mpsB)
            Cstorage{j} = updateCright(Cright, A, [], B);
        end
    end

    %  Convergence check 
    if std(Evalues) / abs(mean(Evalues)) < precision
        % Absorb final gauge into site 1
        mps{1} = contracttensors(mps{1}, 3, 2, U, 2, 1);
        mps{1} = permute(mps{1}, [1, 3, 2]);
        break;
    end
end


% ==========================================================================
function [A, E] = minimizeE_onesite(hsetj, Hleft, Hright, P)
% Solve single-site eigenvalue problem (smallest eigenvalue).

DAl = size(Hleft{1}, 1);
DAr = size(Hright{1}, 1);
d   = size(hsetj{1}, 1);
M   = size(hsetj, 1);

Heff = 0;
for m = 1:M
    % Build effective Hamiltonian matrix in the local basis
    Heffm = contracttensors(Hleft{m},  3, 2, Hright{m}, 3, 2);
    Heffm = contracttensors(Heffm,     5, 5, hsetj{m},  3, 3);
    Heffm = permute(Heffm, [1, 3, 5, 2, 4, 6]);
    Heffm = reshape(Heffm, [DAl*DAr*d, DAl*DAr*d]);
    Heff  = Heff + Heffm;
end

% Project onto orthogonal subspace if targeting excited states
if ~isempty(P), Heff = P' * Heff * P; end

% Find smallest-real eigenvalue and eigenvector
options.disp = 0;
[A, E] = eigs(Heff, 1, 'sr', options);

% Unproject and reshape back to 3-index MPS tensor
if ~isempty(P), A = P * A; end
A = reshape(A, [DAl, DAr, d]);


% ==========================================================================
function [P] = calcprojector_onesite(B, Cleft, Cright)
% Build the projector onto the complement of mpsB.

y = contracttensors(Cleft, 3, 3, B, 3, 1);
y = contracttensors(y, 4, [2, 3], Cright, 3, [2, 3]);
y = permute(y, [1, 3, 2]);
y = reshape(y, [prod(size(y)), 1]);
% Orthonormal basis for the complement: drop the first column (= y direction)
Q = orth([y, eye(size(y, 1))]);
P = Q(:, 2:end);


% ==========================================================================
function [Hstorage] = initHstorage(mps, hset, d)
% Initialize right-environment storage for the Hamiltonian.

[M, N]    = size(hset);
Hstorage  = cell(M, N + 1);
for m = 1:M
    Hstorage{m, 1}   = 1;   % trivial left boundary
    Hstorage{m, N+1} = 1;   % trivial right boundary
end
for j = N:-1:2
    for m = 1:M
        h = reshape(hset{m, j}, [1, 1, d, d]);
        Hstorage{m, j} = updateCright(Hstorage{m, j+1}, mps{j}, h, mps{j});
    end
end
