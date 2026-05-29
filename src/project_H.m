function [x, y] = project_H(Ecut, Emax, eta, step, sn)
% PROJECT_H  Compute the spectral function <O>(E) over an energy window.
%
%   [x, y] = project_H(Ecut, Emax, eta, step, sn)
%
%   Traces out the ratio DOS_O(E) / DOS(E) — i.e. the thermal/microcanonical
%   expectation value of an observable O as a function of energy E — over
%   the spectral range [Ecut, Emax].
%
%   TWO-PHASE APPROACH
%   The spectrum is divided into a "bulk" region (handled by standard KPM
%   via dos.m) and a "fine" region near the spectral edges or sharp features
%   (handled by the higher-resolution dos_precise.m). The transition between
%   phases is triggered when the DOS drops below a fraction eta of its
%   value at the first fine-region boundary.
%
%   Phase 1 (bulk):  use dos.m while DOS > eta * DOS(es)
%   Phase 2 (edges): use dos_precise.m with a sliding reference window
%
%   NOTE
%     This function uses several variables (N, M, mu, muo, nu, del, R) from
%     the caller's workspace via implicit scoping. In a refactored version
%     these should be passed as explicit arguments.
%
%   INPUTS
%     Ecut - Starting energy for the scan.
%     Emax - Upper energy bound (scan stops when E > Emax - 1).
%     eta  - DOS threshold fraction for switching to fine-resolution mode.
%            Typical value: 0.01–0.1.
%     step - Energy step size for the scan grid.
%     sn   - Window size for the sliding reference in Phase 2.
%
%   OUTPUTS
%     x - Energy values where the ratio was evaluated.
%     y - Corresponding values of DOS_O(E) / DOS(E) = <O>_microcanonical(E).
%
%   SEE ALSO: dos, dos_precise, isingtr

e  = Ecut;
es = Ecut;

% First point
x = Ecut;
y = dos(e, N, M, muo, nu, del) / dos(e, N, M, mu, nu, del);
w = dos(es, N, M, mu, nu, del);   % reference DOS value at current window start

% ===== Phase 1: bulk region — use standard KPM (dos.m) ===================
while 1
    e = e + step;
    p = dos(e, N, M, muo, nu, del);   % numerator DOS (with observable)
    q = dos(e, N, M, mu,  nu, del);   % denominator DOS (total)

    if q >= eta * w
        % DOS is still above threshold: record this point
        x = [x, e];
        y = [y, p / q];
    else
        % DOS dropped below threshold: switch to fine-resolution mode
        es = e - step;
        e  = es;
        break;
    end
end

% ===== Phase 2: edge region — use dos_precise.m with sliding window =======
% The sliding window [es, es+sn] is advanced as the DOS continues to decay.
while 1
    % Reference DOS using precise window starting at es
    w = dos_precise(es, es + sn, M, R, mu, N, nu, del);

    while 1
        e = e + step;
        p = dos_precise(e, es + sn, M, R, muo, N, nu, del);
        q = dos_precise(e, es + sn, M, R, mu,  N, nu, del);

        if q >= eta * w
            x = [x, e];
            y = [y, p / q];
        else
            % Advance the reference window
            es = e - step;
            e  = es;
            break;
        end
    end

    % Stop when we have covered the desired energy range
    if es > Emax - 1
        break;
    end
end
