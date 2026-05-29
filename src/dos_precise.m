function Dos = dos_precise(E, Eth, M, R, mu, N, nu, del)
% DOS_PRECISE  Spectral-resolved DOS with two-variable Chebyshev expansion.
%
%   Dos = dos_precise(E, Eth, M, R, mu, N, nu, del)
%
%   Computes a "local" or "energy-windowed" density of states by combining
%   a standard KPM expansion in the energy variable E with a separate
%   spectral resolution expansion at threshold energy Eth. This implements
%   the two-variable Chebyshev expansion:
%
%       D(E, Eth) = sum_{n,m} g_n^M * T_n(Etilde) * g_m^R * alpha_m(Eth_tilde)
%                    * [(mu_{n+m} + mu_{|n-m|})/2] / (2^N * pi * nu * sqrt(...))
%
%   The function alpha_m(x) = sin(m * acos(x)) / (m*pi) for m>0, and
%   alpha_0(x) = acos(x)/pi, which are the Chebyshev coefficients of the
%   step function theta(x - Eth).
%
%   This is useful for computing spectral quantities that require a sharp
%   energy cutoff, such as the projected DOS below a given energy.
%
%   INPUTS
%     E   - Energy values at which to evaluate the spectral function.
%     Eth - Threshold energy (boundary of the spectral window).
%     M   - Number of Chebyshev moments in the energy direction.
%     R   - Number of Chebyshev moments in the resolution direction.
%     mu  - Vector of Chebyshev moments (length >= M+R).
%     N   - System size (number of sites).
%     nu  - Spectral half-width (rescaling parameter).
%     del - Spectral shift (rescaling parameter).
%
%   OUTPUT
%     Dos - Vector of spectral function values at each energy in E.
%
%   SEE ALSO: dos, ker_jackson (local), alpha (local)

L   = length(E);
Dos = [];

for j = 1:L
    % Rescale both energies into the Chebyshev domain [-1, 1]
    Etilde = E(j) / nu + del;
    tEth   = Eth / nu + del;

    % Pre-compute Chebyshev polynomial values T_n(Etilde) for n=0..M-1
    TE = [];
    for i = 1:M
        t  = chebyshevT(i - 1, Etilde);
        TE = [TE, t];
    end

    % Prefactor: 1 / (2^N * pi * nu * sqrt(1 - Etilde^2))
    % (same as the standard KPM normalization)

    % ----- n=0, m=0 term -----
    e = ker_jackson(0, M) * ker_jackson(0, R) * alpha(0, tEth) * mu(1) ...
        / (2^N * pi * nu * sqrt(1 - Etilde^2));

    % ----- n>0, m=0 terms: standard KPM sum for fixed Eth window -----
    for n = 1:(M - 1)
        e = e + 2 * ker_jackson(0, R) * alpha(0, tEth) * ker_jackson(n, M) ...
              * TE(n+1) * mu(n+1) / (2^N * pi * nu * sqrt(1 - Etilde^2));
    end

    % ----- n=0, m>0 terms: resolution expansion at Eth -----
    for m = 1:(R - 1)
        e = e + 2 * ker_jackson(0, M) * alpha(m, tEth) * ker_jackson(m, R) ...
              * mu(m+1) / (2^N * pi * nu * sqrt(1 - Etilde^2));
    end

    % ----- n>0, m>0 cross terms -----
    % Uses the Chebyshev product identity:
    %   T_n(x)*T_m(x) = (T_{n+m}(x) + T_{|n-m|}(x)) / 2
    % leading to the combination (mu_{n+m} + mu_{|n-m|}) / 2
    for n = 1:(M - 1)
        for m = 1:(R - 1)
            e = e + 4 * ker_jackson(n, M) * TE(n+1) * ker_jackson(m, R) ...
                  * alpha(m, tEth) * ((mu(n+m+1) + mu(abs(n-m)+1)) / 2) ...
                  / (2^N * pi * nu * sqrt(1 - Etilde^2));
        end
    end

    Dos = [Dos, e];
end


% =========================================================================
function [gamma] = ker_jackson(n, M)
% KER_JACKSON  Jackson kernel damping coefficient for order n out of M.
%
%   Suppresses Gibbs oscillations by smoothly tapering high-order moments.
%   g_0 = 1 exactly; for n>0:
%       g_n = [(M-n+1)*cos(pi*n/(M+1)) + sin(pi*n/(M+1))*cot(pi/(M+1))] / (M+1)
if n == 0
    gamma = 1;
else
    gamma = ((M - n + 1) * cos(pi*n/(M+1)) + ...
              sin(pi*n/(M+1)) * cot(pi/(M+1))) / (M + 1);
end


% =========================================================================
function [al] = alpha(n, tEth)
% ALPHA  Chebyshev expansion coefficient of the Heaviside step function.
%
%   alpha_0(x) = acos(x) / pi        (fraction of spectrum below x)
%   alpha_n(x) = sin(n*acos(x))/(n*pi)  for n > 0
%
%   These are the expansion coefficients such that
%       theta(x - Eth) = sum_n alpha_n(Eth) * T_n(x)
%   which allows a smooth spectral cutoff in the two-variable expansion.
if n == 0
    al = acos(tEth) / pi;
else
    al = sin(n * acos(tEth)) / (n * pi);
end
