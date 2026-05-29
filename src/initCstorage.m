function [Cstorage] = initCstorage(mpsB, mpoX, mpsA, N)
% INITCSTORAGE  Initialize the right-to-left environment tensor array.
%
%   Cstorage = initCstorage(mpsB, mpoX, mpsA, N)
%
%   Precomputes and stores all right-environment (C) tensors needed for a
%   left-to-right optimization sweep. The environment at site i represents
%   the partial contraction of the network from site i+1 to N:
%
%       Cstorage{i} encodes  <B_{i..N}| [X_{i..N}] |A_{i..N}>
%
%   where B is the bra MPS, A is the ket MPS, and X is an optional MPO.
%   If X is empty, the overlap <B|A> is computed instead.
%
%   This initialization scans from right to left (site N down to site 2),
%   building up the right environments one site at a time using updateCright.
%
%   STORAGE LAYOUT
%     Cstorage{1}   = 1    (trivial left boundary)
%     Cstorage{j}   = right environment from site j to N
%     Cstorage{N+1} = 1    (trivial right boundary)
%
%   INPUTS
%     mpsB  - Cell array {1,N}: bra MPS tensors (typically the MPS being
%             optimized, e.g. the truncated approximation).
%     mpoX  - Cell array {1,N}: MPO tensors. Pass [] to compute overlaps
%             without an operator (identity is used automatically).
%     mpsA  - Cell array {1,N}: ket MPS tensors (e.g. the original MPO-MPS
%             being approximated).
%     N     - Number of lattice sites.
%
%   OUTPUT
%     Cstorage - Cell array {1,N+1} of environment tensors. Each entry has
%                shape [Dl_B, Dl_X, Dl_A] (3-index tensor). Boundary entries
%                are scalars equal to 1.
%
%   SEE ALSO: updateCright, updateCleft, reduceD, minimizeE

Cstorage      = cell(1, N + 1);
Cstorage{1}   = 1;     % trivial left boundary (no sites to the left)
Cstorage{N+1} = 1;     % trivial right boundary (no sites to the right)

% Build right environments from site N down to site 2
for i = N:-1:2
    if isempty(mpoX)
        X = [];             % updateCright handles the identity case internally
    else
        X = mpoX{i};
    end
    Cstorage{i} = updateCright(Cstorage{i+1}, mpsB{i}, X, mpsA{i});
end
