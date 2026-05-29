function [Cstorage] = initCstorage(mpsB, mpoX, mpsA, N)
% Initialize the right-to-left environment tensor array.
%
%   INPUTS
%     mpsB  - Cell array {1,N}: bra MPS tensors.
%     mpoX  - Cell array {1,N}: MPO tensors.
%     mpsA  - Cell array {1,N}: ket MPS tensors.
%     N     - Number of lattice sites.
%
%   OUTPUT
%     Cstorage - Cell array {1,N+1} of environment tensors.
%


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
