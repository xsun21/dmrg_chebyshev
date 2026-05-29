function [Dos] = dos(E, N, M, mu, nu, del)
% DOS  Density of states via the Kernel Polynomial Method (KPM).
%
%   Dos = dos(E, N, M, mu, nu, del)
%
%   Reconstructs the many-body density of states from Chebyshev moments
%   using the Jackson kernel to suppress Gibbs oscillations.

L   = length(E);
Dos = [];

for i = 1:L
    Etilde = E(i) / nu + del;

    gamma = 1;
    e = gamma * mu(1) / (2^N * pi * nu * sqrt(1 - Etilde^2));

    for n = 1:(M - 1)
        gamma = ((M - n + 1) * cos(pi*n/(M+1)) + ...
                  sin(pi*n/(M+1)) * cot(pi/(M+1))) / (M + 1);
        TE = chebyshevT_eval(n, Etilde);
        e  = e + 2 * gamma * mu(n+1) * TE / (2^N * pi * nu * sqrt(1 - Etilde^2));
    end

    Dos = [Dos, e];
end
