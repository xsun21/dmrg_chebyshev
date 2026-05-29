function [X, numindX] = contracttensors(X, numindX, indX, Y, numindY, indY)
% CONTRACTTENSORS  Contract two tensors along specified index pairs.
%
%   [X, numindX] = contracttensors(X, numindX, indX, Y, numindY, indY)
%
%   This is the core tensor contraction primitive used throughout the MPS/MPO
%   framework. It generalizes matrix multiplication to arbitrary-order tensors.
%   Conceptually it computes:
%
%       Z_{i1,...,ik, j1,...,jl} = sum_{c1,...,cm} X_{...,c1,...,cm,...} *
%                                                   Y_{c1,...,cm,...}
%
%   where indX and indY specify which indices are summed (contracted) over.
%
%   The implementation works by:
%     1. Permuting contracted indices to the right of X and left of Y.
%     2. Reshaping both tensors into matrices.
%     3. Performing standard matrix multiplication.
%     4. Reshaping the result back into a tensor.
%
%   INPUTS
%     X       - First input tensor (any order array).
%     numindX - Number of indices (logical order) of X. Must be >= ndims(X).
%               Allows treating lower-dimensional arrays as higher-order tensors
%               with trailing singleton dimensions.
%     indX    - Vector of indices in X to be contracted (summed over).
%     Y       - Second input tensor.
%     numindY - Number of indices of Y.
%     indY    - Vector of indices in Y to be contracted. Must match indX in
%               total size: prod(size_X(indX)) == prod(size_Y(indY)).
%
%   OUTPUTS
%     X       - Result tensor containing the uncontracted (free) indices of
%               X followed by the free indices of Y.
%     numindX - Number of indices of the result.
%
%   EXAMPLE
%     % Matrix-vector product: C = A * b  (contract index 2 of A with index 1 of b)
%     A = rand(3,4); b = rand(4,1);
%     [C, ~] = contracttensors(A, 2, 2, b, 2, 1);  % C is 3x1

% Read out sizes, padding to numind in case of trailing singleton dimensions
Xsize = ones(1, numindX);  Xsize(1:length(size(X))) = size(X);
Ysize = ones(1, numindY);  Ysize(1:length(size(Y))) = size(Y);

% Free (uncontracted) index lists
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

% ---- Special case: Y has no free indices (result is a scalar or vector) ----
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

% ---- General case: both tensors have free indices -------------------------
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
