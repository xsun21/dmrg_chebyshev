function Dos = dos_precise(E, Eth, M, R, mu, N, nu, del)
% Spectral-resolved DOS with two-variable Chebyshev expansion.
%
%   INPUTS
%     E   - Energy values at which to evaluate the spectral function.
%     Eth - Threshold energy (boundary of the spectral window).
%     M   - Number of Chebyshev moments in the energy direction.
%     R   - Number of Chebyshev moments in the resolution direction.
%     mu  - Vector of Chebyshev moments.
%     N   - System size.
%     nu  - Spectral half-width (rescaling parameter).
%     del - Spectral shift (rescaling parameter).
%
%   OUTPUT
%     Dos - Vector of spectral function values at each energy in E.
%

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
% Jackson kernel damping coefficient for order n out of M.
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
% Chebyshev expansion coefficient of the Heaviside step function.
if n == 0
    al = acos(tEth) / pi;
else
    al = sin(n * acos(tEth)) / (n * pi);
end
