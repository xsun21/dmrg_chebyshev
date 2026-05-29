% ISINGTR  Full density-of-states calculation for the Ising model via KPM + MPS.
%
%   Computes the density of states rho(E) and the energy-resolved expectation
%   value <sigma^z_{N/2}>(E) for the 1D mixed-field Ising Hamiltonian using
%   the Chebyshev expansion of the MPO T_n(H_tilde) up to order M.
%
%   WORKFLOW
%     1. Rescale H -> H_tilde = (H - b*I) / nu  so spectrum lies in [-1, 1].
%     2. Build the MPO representation of H_tilde.
%     3. Call ChebyshevH to compute Chebyshev moments mu_n = Tr[T_n(H_tilde)]
%        and (optionally) muo_n = Tr[O * T_n(H_tilde)] for observable O.
%     4. Reconstruct rho(E) via the KPM formula in dos.m.
%     5. Plot <O>(E) = DOS_O(E) / DOS(E) vs E.
%
%   SPECTRAL RESCALING
%     The Chebyshev recursion requires the Hamiltonian's eigenvalues to lie
%     strictly within (-1, 1). Given spectral bounds Emin and Emax (from DMRG):
%
%       ep   = 1e-8                   (small margin to avoid endpoints)
%       nu   = (Emax-Emin) / (2-ep)   (half-width)
%       del  = -((2-ep)*(Emax+Emin)) / (2*(Emax-Emin))   (shift)
%
%     so H_tilde = H/nu + del*I has eigenvalues in [-(1-ep/2), 1-ep/2].
%
%   MPO FORMAT
%     H_tilde is encoded as an MPO in the standard "W-matrix" form.
%     For the Ising chain the W-matrix at bulk sites is:
%
%         W = [ id     0     0  ]
%             [ a      0     0  ]
%             [ b     sz    id  ]
%
%     where a = J*sz/nu, b = (g*sx + h*sz)/nu + del*id/N.
%     Left/right boundary conditions terminate the virtual indices.
%
%   OUTPUT FILES
%     trH.txt   - Chebyshev moments mu_n = Tr[T_n(H_tilde)]
%     trOH.txt  - Chebyshev moments muo_n = Tr[sigma^z * T_n(H_tilde)]
%
%   SEE ALSO: ChebyshevH, dos, ising, minimizeE

%% --- System parameters --------------------------------------------------
N  = 40;
d  = 2;      % spin-1/2 local dimension
dw = 3;      % MPO bond dimension for H_tilde (3-state W-matrix)
M  = 100;    % number of Chebyshev moments (resolution ~ pi/M in rescaled units)

%% --- Convergence and truncation -----------------------------------------
precision = 1e-8;    % variational compression convergence
trunD     = 200;     % maximum MPO bond dimension before truncation

%% --- Spectral bounds (obtained from ising.m DMRG) -----------------------
% These values should be updated if N, J, g, h are changed.
Emin = -52.903434036554550;
Emax =  68.196022812764740;

ep  = 1e-8;
nu  = (Emax - Emin) / (2 - ep);
del = -((2 - ep) * (Emax + Emin)) / (2 * (Emax - Emin));

%% --- Pauli matrices and Hamiltonian couplings ---------------------------
J  = 1;
g  = -1.05;
h  = 0.5;

sx = [0 1; 1 0];
sz = [1 0; 0 -1];
id = [1 0; 0 1];
ze = [0 0; 0 0];   % zero matrix (placeholder)

% Rescaled single-site operators
a = J * sz / nu;                         % ZZ coupling contribution
b = g * sx / nu + h * sz / nu + del * id / N;   % on-site terms (incl. shift)

%% --- Construct MPO for H_tilde ------------------------------------------
% The W-matrix representation of H_tilde has virtual bond dimension dw=3:
%   virtual index 1 = "identity channel" (right boundary)
%   virtual index 2 = "sz channel" (ZZ bond carrying)
%   virtual index 3 = "on-site channel" (left boundary accumulates terms)
%
% Site tensors have shape H{i}(Dl, Dr, d, d).

H = cell(1, N);

% Left boundary site (site 1): generates terms, no left virtual index needed
H{1}(1, 1, :, :) = b;    % on-site term
H{1}(1, 2, :, :) = sz;   % start of ZZ bond: left sz
H{1}(1, 3, :, :) = id;   % pass-through identity

% Right boundary site (site N): terminates all channels
H{N}(1, 1, :, :) = id;   % identity closes the pass-through
H{N}(2, 1, :, :) = a;    % right sz closes the ZZ bond: a = J*sz/nu
H{N}(3, 1, :, :) = b;    % on-site term

% Bulk sites (2 <= i <= N-1): propagate all three channels
for i = 2:(N - 1)
    H{i}(1, 1, :, :) = id;   % identity pass-through
    H{i}(2, 1, :, :) = a;    % ZZ bond termination
    H{i}(3, 1, :, :) = b;    % on-site accumulation
    H{i}(3, 2, :, :) = sz;   % start new ZZ bond
    H{i}(3, 3, :, :) = id;   % pass-through
end

%% --- Construct MPO for observable: sigma^z at site N/2 ------------------
% The observable is a single-site operator, so its MPO is diagonal with
% identity everywhere except at site N/2 where it equals sz.

O = cell(1, N);
for i = 1:N
    O{i}(1, 1, :, :) = id;   % identity at all sites
end
O{N/2}(1, 1, :, :) = sz;     % observable at center site

%% --- Chebyshev expansion ------------------------------------------------
% Computes:
%   mu(n)  = Tr[ T_{n-1}(H_tilde) ]          (partition function moments)
%   muo(n) = Tr[ O * T_{n-1}(H_tilde) ]      (observable moments)
% for n = 1, 2, ..., M+1 (orders 0 through M).

[mu, muo, err] = ChebyshevH(H, O, 200, trunD, precision);

%% --- Save moments to disk -----------------------------------------------
dlmwrite('trH.txt',  mu,  'precision', '%.10f')
dlmwrite('trOH.txt', muo, 'precision', '%.10f')

%% --- Compute and plot DOS ratio <sigma^z>(E) ----------------------------
% Evaluate both DOS curves on a coarse energy grid
x  = -20:2:70;
y1 = dos(x, N, M, muo, nu, del);   % DOS weighted by <sigma^z>
y2 = dos(x, N, M, mu,  nu, del);   % total DOS

% Ratio = microcanonical expectation value of sigma^z at energy E
z = length(x);
l = zeros(1, z);
for i = 1:z
    l(i) = y1(i) / y2(i);
end

figure;
plot(x, l);
xlabel('Energy E');
ylabel('\langle\sigma^z_{N/2}\rangle(E)');
title('Microcanonical expectation value vs energy (Ising model, N=40)');
grid on;
