function [x, y] = project_H(Ecut, Emax, eta, step, sn)
% PROJECT_H  Compute the spectral function <O>(E) over an energy window.
%
%   [x, y] = project_H(Ecut, Emax, eta, step, sn)
%
%   Traces out the ratio DOS_O(E) / DOS(E) — i.e. the thermal/microcanonical
%   expectation value of an observable O as a function of energy E — over
%   the spectral range [Ecut, Emax].

e  = Ecut;
es = Ecut;

% First point
x = Ecut;
y = dos(e, N, M, muo, nu, del) / dos(e, N, M, mu, nu, del);
w = dos(es, N, M, mu, nu, del);   % reference DOS value at current window start

% ===== Phase 1: bulk region — use standard KPM (dos.m) =====
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

% ===== Phase 2: edge region - use dos_precise.m with sliding window =======
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
