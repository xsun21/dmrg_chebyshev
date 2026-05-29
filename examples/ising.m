% ISING  DMRG ground and maximum eigenstate for the transverse-field Ising model.
%
%   Computes the ground-state energy E0 and maximum energy E1 of the
%   1D mixed-field Ising Hamiltonian:
%
%       H = J * sum_{i} sz_i * sz_{i+1}
%         + g * sum_{i} sx_i
%         + h * sum_{i} sz_i
%
%   PARAMETERS
%     N   = 40   lattice sites
%     D   = 14   MPS bond dimension (increase for better accuracy)
%     J   = 1    Ising coupling
%     g   = -1.05 transverse field
%     h   = 0.5   longitudinal field
%
%   OUTPUT
%     Prints: D=14, E0=<ground state energy>, E1=<max energy>
%

%% --- Parameters ---
N         = 40;
D         = 14;
J         = 1;
g         = -1.05;
h         = 0.5;
precision = 1e-10;

%% --- Pauli matrices and building blocks ---
sx = [0, 1; 1, 0];
sy = [0, -1i; 1i, 0];
sz = [1, 0; 0, -1];
id = eye(2);

%% --- Construct Hamiltonian in hset format ---
% H has three types of terms:
%   (1) N-1 two-site terms: J * sz_j * sz_{j+1}   (j = 1..N-1)
%   (2) N   one-site terms: g * sx_j              (j = 1..N)
%   (3) N   one-site terms: h * sz_j              (j = 1..N)
% Total: M = (N-1) + N + N = 3N-1 terms

M    = N - 1 + 2*N;
hset = cell(M, N);

% Initialize all entries to identity (so that unused sites are trivial)
for m = 1:M
    for j = 1:N
        hset{m, j} = id;
    end
end

% (1) Two-site ZZ coupling: hset{j, j} = J*sz, hset{j, j+1} = sz
for j = 1:(N - 1)
    hset{j, j}   = J * sz;
    hset{j, j+1} = sz;
end

% (2) Transverse field: hset{j+N-1, j} = g*sx
for j = N:(2*N - 1)
    hset{j, j - N + 1} = g * sx;
end

% (3) Longitudinal field: hset{j+2N-1, j} = h*sz
for j = (2*N):(3*N - 1)
    hset{j, j - 2*N + 1} = h * sz;
end

%% --- DMRG optimization ---
randn('state', 0);   % fix random seed for reproducibility

% Ground state (minimum energy)
[E0, mps0] = minimizeE(hset, D, precision, []);

% Maximum energy eigenstate (needed for spectral bounds in isingtr.m)
[E1, mps1] = maximizeE(hset, D, precision, []);

fprintf('D=%g, E0 = %.6g, E1 = %.6g\n', D, E0, E1);
