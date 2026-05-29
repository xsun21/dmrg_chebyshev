function [mu, muo, err, TnH] = ChebyshevH(H, mpo, n, trunD, precision)
% Build the Chebyshev MPO expansion T_n(H) and compute moments.
%
%   INPUTS
%     H         - Cell array {1,N} representing the MPO of the (rescaled) Hamiltonian
%     mpo       - Optional MPO for an observable O. 
%     n         - Maximum Chebyshev order to compute (inclusive).
%     trunD     - Bond dimension threshold.
%     precision - Convergence tolerance passed to reduceD (e.g. 1e-8).
%
%   OUTPUTS
%     mu        - Row vector of length n+1 with moments mu(k) = Tr[T_{k-1}(H)].
%     muo       - Row vector of moments Tr[O * T_{k-1}(H)] (empty if mpo=[]).
%     err       - Accumulated truncation error (sum of relative HS-norm errors).
%     TnH       - Cell MPO representing T_n(H) at the final order.

N  = size(H, 2);          % number of lattice sites
dw = size(H{1}, 2);       % MPO bond dimension of H (right)
d  = size(H{1}, 3);       % local physical dimension (d=2 for spin-1/2)

% Allocate cell arrays for the current (T1) and previous (T0) Chebyshev MPOs
TnH    = cell(1, N);
trunTnH = cell(1, N);
T0     = cell(1, N);   % T_{n-1}(H)  — starts as identity
T1     = cell(1, N);   % T_n(H)      — starts as H

id = [1 0; 0 1];   % 2x2 identity
ze = [0 0; 0 0];   % 2x2 zero matrix (used as placeholder)

mu  = [];   % moments for Tr[T_n(H)]
muo = [];   % moments for Tr[O * T_n(H)]
err = 0;    % cumulative truncation error

% Initialize T0 = Identity MPO
for i = 1:N
    T0{i}(1, 1, :, :) = id;
end

% T1 = H (the Hamiltonian MPO is already the degree-1 Chebyshev polynomial)
T1 = H;

% Compute the zeroth and first moments
mu0 = trChebyshevH(T0, []);   % Tr[T_0(H)] = Tr[I] = 2^N (before normalization)
mu1 = trChebyshevH(T1, []);   % Tr[T_1(H)] = Tr[H]

if ~isempty(mpo)
    muo0 = trChebyshevH(T0, mpo);
    muo1 = trChebyshevH(T1, mpo);
end

% Base cases
if n == 0
    TnH = T0;
    mu  = [mu0];
    if ~isempty(mpo), muo = [muo0]; end
end

if n == 1
    TnH = T1;
    mu  = [mu0, mu1];
    if ~isempty(mpo), muo = [muo0, muo1]; end
end

% Recursive case: apply T_{j+1} = 2*H*T_j - T_{j-1} for j = 1 ... n-1
if n >= 2
    mu = [mu0, mu1];
    if ~isempty(mpo), muo = [muo0, muo1]; end

    fprintf('0,%g\n', mu0);
    fprintf('1,%g\n', mu1);

    for j = 1:(n - 1)
        e = 0;

        % --- Step 1: Compute 2*H*T_j(H) as a new MPO ---
        for i = 1:N
            Hl  = size(H{i},  1);   Hr  = size(H{i},  2);
            lT1 = size(T1{i}, 1);   rT1 = size(T1{i}, 2);

            % Contract over shared physical index (index 4 of H with index 3 of T1)
            % Result has combined bond dims [Hl*lT1, Hr*rT1, d, d]
            TnH{i} = contracttensors(H{i}, 4, 4, T1{i}, 4, 3);
            TnH{i} = permute(TnH{i}, [1, 4, 2, 5, 3, 6]);
            TnH{i} = reshape(TnH{i}, [Hl*lT1, Hr*rT1, d, d]);
        end

        % Multiply the leftmost site by 2 (implements the factor of 2 in 2H*T_j)
        TnH{1} = 2 * TnH{1};

        % --- Step 2: Subtract T_{j-1}(H) = T0 ---
        % MPO subtraction is done by block-diagonal construction in the virtual (bond) indices
        %
        % For the leftmost site: append T0 columns to the right bond.
        T0{1} = -T0{1};
        THr  = size(TnH{1}, 2);   rT0 = size(T0{1}, 2);
        TnH{1}(1, (THr+1):(THr+rT0), :, :) = T0{1};

        % For the rightmost site: append T0 rows to the left bond.
        THl  = size(TnH{N}, 1);   lT0 = size(T0{N}, 1);
        TnH{N}((THl+1):(THl+lT0), 1, :, :) = T0{N};

        % For bulk sites: block-diagonal extension in both bond directions.
        for i = 2:(N - 1)
            THl = size(TnH{i}, 1);   THr = size(TnH{i}, 2);
            lT0 = size(T0{i}, 1);    rT0 = size(T0{i}, 2);
            TnH{i}((THl+1):(THl+lT0), (THr+1):(THr+rT0), :, :) = T0{i};
        end

        % --- Step 3: Truncate bond dimension if necessary ---
        THr = size(TnH{N/2}, 2);   % check bond dim at the center cut
        if THr > trunD
            % Reshape to 3-index tensors (merge physical dims) for reduceD
            for i = 1:N
                THl = size(TnH{i}, 1);   THr = size(TnH{i}, 2);
                TnH{i} = reshape(TnH{i}, [THl, THr, d*d]);
            end

            % Variationally compress TnH to bond dimension trunD
            trunTnH = reduceD(TnH, [], trunD, precision);

            % Restore 4-index form [Dl, Dr, d, d] after compression
            for i = 1:N
                THl  = size(TnH{i},    1);   THr  = size(TnH{i},    2);
                tTHl = size(trunTnH{i}, 1);  tTHr = size(trunTnH{i}, 2);
                TnH{i}    = reshape(TnH{i},    [THl,  THr,  d, d]);
                trunTnH{i} = reshape(trunTnH{i}, [tTHl, tTHr, d, d]);
            end

            % Accumulate relative truncation error in the HS-norm sense:
            %   err += 1 - 2*<trun,full> / <full,full>
            % (ideally 0 if truncation is exact)
            err = err + 1 + (HSnorm(trunTnH, trunTnH) ...
                           - HSnorm(trunTnH, TnH) ...
                           - HSnorm(TnH, trunTnH)) / HSnorm(TnH, TnH);
            TnH = trunTnH;
        end

        % --- Step 4: Compute and store moments ---
        e   = trChebyshevH(TnH, []);
        mu  = [mu, e];
        fprintf('%g,%g\n', j+1, e);

        if ~isempty(mpo)
            e   = trChebyshevH(TnH, mpo);
            muo = [muo, e];
        end

        % Advance recurrence: T_{j-1} <- T_j,  T_j <- T_{j+1}
        T0 = T1;
        T1 = TnH;
    end
end
