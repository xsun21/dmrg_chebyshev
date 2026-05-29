function [mps] = prepare(mps)
% Bring an MPS into right-canonical form.
%   INPUT / OUTPUT
%     mps - Cell array {1,N} of MPS tensors, each of shape [Dl, Dr, d].


N = length(mps);
for i = N:-1:2
    % Right-canonical decompose site i:
    [mps{i}, U] = prepare_onesite(mps{i}, 'rl');

    % Absorb gauge matrix U into the left neighbor:
    %   mps{i-1}[:, :, s] <- mps{i-1}[:, :, s] * U
    % using index 2 of mps{i-1} (right bond) contracted with index 1 of U
    mps{i-1} = contracttensors(mps{i-1}, 3, 2, U, 2, 1);
    mps{i-1} = permute(mps{i-1}, [1, 3, 2]);   % restore [Dl, Dr, d] layout
end
