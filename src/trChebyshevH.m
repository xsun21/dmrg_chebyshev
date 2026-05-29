function [tr] = trChebyshevH(TnH, mpo)
% TRCHEBYSHEVH  Compute Tr[T_n(H)] or Tr[O * T_n(H)] from MPO representations.
%
%   tr = trChebyshevH(TnH, mpo)
%
%   Evaluates the many-body trace of a Chebyshev MPO (or its product with
%   an observable MPO) by exploiting the MPO structure:
%
%       Tr[T_n(H)] = product_{i=1}^{N} tr_i[ TnH{i} ]
%
%   where tr_i denotes the trace over the physical indices at site i.
%   This factorizes because the trace of a tensor product equals the product
%   of traces: Tr[A ⊗ B] = Tr[A] * Tr[B].
%
%   If an observable MPO O is provided, the code first computes the product
%   MPO U = O * TnH (by contracting physical indices site by site), then
%   traces U in the same way.
%
%   INPUTS
%     TnH - Cell array {1,N} of Chebyshev MPO tensors, each of shape
%           [Dl, Dr, d, d]. Typically the output of ChebyshevH.
%     mpo - Optional observable MPO (cell array {1,N} of shape [Dl, Dr, d, d]).
%           Pass [] to compute Tr[TnH] directly.
%
%   OUTPUT
%     tr  - Scalar trace value. This is the raw (unnormalized) trace, equal
%           to 2^N * mu_n where mu_n is the n-th Chebyshev moment.
%
%   SEE ALSO: ChebyshevH, dos, dos_precise

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
% Because trs{i} are matrices (virtual bond indices), the full trace is
% their matrix product contracted to a scalar (exploiting open boundaries):
%   tr = trs{1} * trs{2} * ... * trs{N}
% which collapses to a scalar because trs{1} has shape [1, Dr] and trs{N} has [Dl, 1].
tr = 1;
for i = 1:N
    tr = tr * trs{i};
end
