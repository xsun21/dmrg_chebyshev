function [mps] = createrandommps(N, D, d)
% Initialize a random Matrix Product State (MPS).
%
%   INPUTS
%     N - Number of lattice sites.
%     D - Virtual bond dimension.
%     d - Physical dimension per site (d=2 for spin-1/2).
%
%   OUTPUT
%     mps - Cell array {1,N} of random MPS tensors.

mps    = cell(1, N);
mps{1} = randn(1, D, d) / sqrt(D);   % left boundary: bond dim 1 on left
mps{N} = randn(D, 1, d) / sqrt(D);   % right boundary: bond dim 1 on right
for i = 2:(N - 1)
    mps{i} = randn(D, D, d) / sqrt(D);   % bulk sites: full bond dimension
end
