function [mps] = createrandommps(N, D, d)
% CREATERANDOMMPS  Initialize a random Matrix Product State (MPS).
%
%   mps = createrandommps(N, D, d)
%
%   Creates a cell array representing a random MPS with open boundary
%   conditions. Each tensor mps{i} has shape [Dl, Dr, d] where:
%     - Dl is the left virtual bond dimension
%     - Dr is the right virtual bond dimension
%     - d  is the physical (local Hilbert space) dimension
%
%   Boundary sites have bond dimension 1 on the open end:
%     - mps{1}:  [1, D, d]
%     - mps{N}:  [D, 1, d]
%     - mps{i}:  [D, D, d]  for 2 <= i <= N-1
%
%   Entries are drawn from a standard normal distribution and divided by
%   sqrt(D) so that the norm of the state stays O(1) for moderate D.
%
%   This state is NOT normalized or in any canonical form. Call prepare()
%   afterward to bring it into right-canonical form before use in DMRG or
%   variational compression.
%
%   INPUTS
%     N - Number of lattice sites.
%     D - Virtual bond dimension (controls entanglement capacity).
%     d - Physical dimension per site (d=2 for spin-1/2).
%
%   OUTPUT
%     mps - Cell array {1,N} of random MPS tensors.
%
%   SEE ALSO: prepare, minimizeE, reduceD

mps    = cell(1, N);
mps{1} = randn(1, D, d) / sqrt(D);   % left boundary: bond dim 1 on left
mps{N} = randn(D, 1, d) / sqrt(D);   % right boundary: bond dim 1 on right
for i = 2:(N - 1)
    mps{i} = randn(D, D, d) / sqrt(D);   % bulk sites: full bond dimension
end
