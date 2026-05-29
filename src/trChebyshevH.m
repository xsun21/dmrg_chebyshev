function [tr] = trChebyshevH(TnH, mpo)
% Compute Tr[T_n(H)] or Tr[O * T_n(H)] from MPO representations.
%
%   INPUTS
%     TnH - Cell array {1,N} of Chebyshev MPO tensors, each of shape [Dl, Dr, d, d].
%     mpo - Optional observable MPO (cell array {1,N} of shape [Dl, Dr, d, d]).
%
%   OUTPUT
%     tr  - Scalar trace value. 
%

N = length(TnH);
d = size(TnH{1}, 3);   % physical dimension

U    = cell(1, N);   % product MPO (only used when mpo is non-empty)
trs  = cell(1, N);   % per-site partial traces

if isempty(mpo)
    % ----- Direct trace of TnH -----
    % For each site i, compute the matrix trs{i}(j,k) = sum_s TnH{i}(j,k,s,s)
    % i.e., trace over the physical indices while keeping the virtual (bond) indices.
    for i = 1:N
        THl = size(TnH{i}, 1);
        THr = size(TnH{i}, 2);
        for j = 1:THl
            for k = 1:THr
                e = 0;
                for m = 1:d
                    e = e + TnH{i}(j, k, m, m);   % sum diagonal physical elements
                end
                trs{i}(j, k) = e;
            end
        end
    end

else
    % ----- Trace of (mpo * TnH) -----
    % First build the product MPO U{i} = mpo{i} ⊗_phys TnH{i}:
    % contract the ket physical index of mpo with the bra physical index of TnH.
    for i = 1:N
        THl = size(TnH{i}, 1);   THr = size(TnH{i}, 2);
        Xl  = size(mpo{i},  1);  Xr  = size(mpo{i},  2);

        % Contract index 4 of TnH (ket physical) with index 3 of mpo (bra physical)
        U{i} = contracttensors(TnH{i}, 4, 4, mpo{i}, 4, 3);
        U{i} = permute(U{i}, [1, 4, 2, 5, 3, 6]);
        % Result shape: [THl*Xl, THr*Xr, d, d]  (combined bond dims)
        U{i} = reshape(U{i}, [THl*Xl, THr*Xr, d, d]);
    end

    % Then trace the product MPO U over physical indices
    for i = 1:N
        Ul = size(U{i}, 1);
        Ur = size(U{i}, 2);
        for j = 1:Ul
            for k = 1:Ur
                e = 0;
                for m = 1:d
                    e = e + U{i}(j, k, m, m);
                end
                trs{i}(j, k) = e;
            end
        end
    end
end

% ----- Final contraction: multiply all site-trace matrices -----
tr = 1;
for i = 1:N
    tr = tr * trs{i};
end
