function [E, mps] = maximizeE(hset, D, precision, mpsB)
% MAXIMIZEE  DMRG maximum eigenstate search (single-site, 1D).
%
%   [E, mps] = maximizeE(hset, D, precision, mpsB)
%
%   Identical in structure to minimizeE, but targets the LARGEST eigenvalue
%   instead of the smallest. This is used to find the upper edge of the
%   spectrum (Emax), which is required for the spectral rescaling in isingtr.
%
%   The only algorithmic difference from minimizeE is that eigs is called
%   with 'la' (largest algebraic) instead of 'sr' (smallest real).
%
%   INPUTS / OUTPUTS: identical to minimizeE. See minimizeE for full docs.
%
%   SEE ALSO: minimizeE, isingtr

[M, N] = size(hset);
d      = size(hset{1, 1}, 1);

mps = createrandommps(N, D, d);
mps = prepare(mps);

Hstorage = initHstorage(mps, hset, d);
if ~isempty(mpsB)
    Cstorage = initCstorage(mps, [], mpsB, N);
end
P = [];

while 1
    Evalues = [];

    % ===== Left sweep =====================================================
    for j = 1:(N - 1)
        if ~isempty(mpsB)
            B      = mpsB{j};
            Cleft  = Cstorage{j};
            Cright = Cstorage{j+1};
            P      = calcprojector_onesite(B, Cleft, Cright);
        end

        Hleft  = Hstorage(:, j);
        Hright = Hstorage(:, j+1);
        hsetj  = hset(:, j);

        % Use 'la' = largest algebraic eigenvalue (vs 'sr' in minimizeE)
        [A, E]  = maximizeE_onesite(hsetj, Hleft, Hright, P);
        [A, U]  = prepare_onesite(A, 'lr');
        mps{j}  = A;
        Evalues = [Evalues, E];

        for m = 1:M
            h = reshape(hset{m, j}, [1, 1, d, d]);
            Hstorage{m, j+1} = updateCleft(Hleft{m}, A, h, A);
        end
        if ~isempty(mpsB)
            Cstorage{j+1} = updateCleft(Cleft, A, [], B);
        end
    end

    % ===== Right sweep ====================================================
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

        [A, E]  = maximizeE_onesite(hsetj, Hleft, Hright, P);
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

    if std(Evalues) / abs(mean(Evalues)) < precision
        mps{1} = contracttensors(mps{1}, 3, 2, U, 2, 1);
        mps{1} = permute(mps{1}, [1, 3, 2]);
        break;
    end
end


% ==========================================================================
function [A, E] = maximizeE_onesite(hsetj, Hleft, Hright, P)
% MAXIMIZEE_ONESITE  Single-site eigenvalue problem (largest eigenvalue).

DAl = size(Hleft{1}, 1);
DAr = size(Hright{1}, 1);
d   = size(hsetj{1}, 1);
M   = size(hsetj, 1);

Heff = 0;
for m = 1:M
    Heffm = contracttensors(Hleft{m},  3, 2, Hright{m}, 3, 2);
    Heffm = contracttensors(Heffm,     5, 5, hsetj{m},  3, 3);
    Heffm = permute(Heffm, [1, 3, 5, 2, 4, 6]);
    Heffm = reshape(Heffm, [DAl*DAr*d, DAl*DAr*d]);
    Heff  = Heff + Heffm;
end

if ~isempty(P), Heff = P' * Heff * P; end

options.disp = 0;
% Key difference vs minimizeE: 'la' = largest algebraic eigenvalue
[A, E] = eigs(Heff, 1, 'la', options);

if ~isempty(P), A = P * A; end
A = reshape(A, [DAl, DAr, d]);


% ==========================================================================
function [P] = calcprojector_onesite(B, Cleft, Cright)
y = contracttensors(Cleft, 3, 3, B, 3, 1);
y = contracttensors(y, 4, [2, 3], Cright, 3, [2, 3]);
y = permute(y, [1, 3, 2]);
y = reshape(y, [prod(size(y)), 1]);
Q = orth([y, eye(size(y, 1))]);
P = Q(:, 2:end);


% ==========================================================================
function [Hstorage] = initHstorage(mps, hset, d)
[M, N]   = size(hset);
Hstorage = cell(M, N + 1);
for m = 1:M
    Hstorage{m, 1}   = 1;
    Hstorage{m, N+1} = 1;
end
for j = N:-1:2
    for m = 1:M
        h = reshape(hset{m, j}, [1, 1, d, d]);
        Hstorage{m, j} = updateCright(Hstorage{m, j+1}, mps{j}, h, mps{j});
    end
end
