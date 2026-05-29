function [mps] = prepare(mps)
% PREPARE  Bring an MPS into right-canonical form.
%
%   mps = prepare(mps)
%
%   Sweeps from site N to site 2, applying right-canonical gauge fixing
%   (via SVD) at each site and absorbing the resulting singular value matrix
%   into the left neighbor. After this call, all tensors mps{i} for i >= 2
%   satisfy the right-canonical condition:
%
%       sum_s mps{i}(:,:,s)' * mps{i}(:,:,s) = I
%
%   The gauge freedom is absorbed entirely into mps{1}, which is left in
%   no particular canonical form (it carries the full norm of the state).
%
%   This initialization is a prerequisite for DMRG sweeps (minimizeE,
%   maximizeE) and variational compression (reduceD), which assume the
%   initial MPS is in right-canonical form before the left-to-right sweep.
%
%   INPUT / OUTPUT
%     mps - Cell array {1,N} of MPS tensors, each of shape [Dl, Dr, d].
%           Modified in-place and returned.
%
%   SEE ALSO: prepare_onesite, minimizeE, reduceD

N = length(mps);
for i = N:-1:2
    % Right-canonical decompose site i: mps{i} -> B (right-canonical) + U (gauge)
    [mps{i}, U] = prepare_onesite(mps{i}, 'rl');

    % Absorb gauge matrix U into the left neighbor:
    %   mps{i-1}[:, :, s] <- mps{i-1}[:, :, s] * U
    % using index 2 of mps{i-1} (right bond) contracted with index 1 of U
    mps{i-1} = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
    mps{i-1} = permute(mps{i-1}, [1, 3, 2]);   % restore [Dl, Dr, d] layout
end
