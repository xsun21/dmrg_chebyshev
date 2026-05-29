function [X, numindX] = contracttensors(X, numindX, indX, Y, numindY, indY)
% Contract two tensors along specified index pairs.
%
%   INPUTS
%     X       - First input tensor (any order array).
%     numindX - Number of indices of X.
%     indX    - Vector of indices in X to be contracted.
%     Y       - Second input tensor.
%     numindY - Number of indices of Y.
%     indY    - Vector of indices in Y to be contracted.
%
%   OUTPUTS
%     X       - Result tensor containing the uncontracted indices.
%     numindX - Number of indices of the result.
%

% Read out sizes, padding to numind in case of trailing singleton dimensions
Xsize = ones(1, numindX);  Xsize(1:length(size(X))) = size(X);
Ysize = ones(1, numindY);  Ysize(1:length(size(Y))) = size(Y);

% Free index lists
indXl = 1:numindX;  indXl(indX) = [];   % free indices of X
indYr = 1:numindY;  indYr(indY) = [];   % free indices of Y

% Sizes of contracted and free index groups
sizeXl = Xsize(indXl);
sizeX  = Xsize(indX);
sizeYr = Ysize(indYr);
sizeY  = Ysize(indY);

% Sanity check: contracted dimensions must match
if prod(sizeX) ~= prod(sizeY)
    error('indX and indY are not of same dimension.');
end

% Special case: Y has no free indices (result is a scalar or vector)
if isempty(indYr)
    if isempty(indXl)
        % Full contraction -> scalar
        X = permute(X, [indX]);
        X = reshape(X, [1, prod(sizeX)]);
        Y = permute(Y, [indY]);
        Y = reshape(Y, [prod(sizeY), 1]);
        X = X * Y;
        Xsize = 1;
        return;
    else
        % X has free indices, Y does not -> result has shape of X's free indices
        X = permute(X, [indXl, indX]);
        X = reshape(X, [prod(sizeXl), prod(sizeX)]);
        Y = permute(Y, [indY]);
        Y = reshape(Y, [prod(sizeY), 1]);
        X = X * Y;
        Xsize = Xsize(indXl);
        X = reshape(X, [Xsize, 1]);
        return
    end
end

% General case: both tensors have free indices 
% Move contracted indices to the right of X
X = permute(X, [indXl, indX]);
X = reshape(X, [prod(sizeXl), prod(sizeX)]);

% Move contracted indices to the left of Y
Y = permute(Y, [indY, indYr]);
Y = reshape(Y, [prod(sizeY), prod(sizeYr)]);

% Matrix multiply: (free-X) x (contract) * (contract) x (free-Y)
X = X * Y;

% Reshape result back to tensor with free indices of X then Y
Xsize    = [Xsize(indXl), Ysize(indYr)];
numindX  = length(Xsize);
X        = reshape(X, [Xsize, 1]);
