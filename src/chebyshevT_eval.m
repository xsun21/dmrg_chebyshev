function T = chebyshevT_eval(n, x)
% CHEBYSHEVT_EVAL  Evaluate Chebyshev polynomial T_n(x) without toolboxes.
%
%   Uses the three-term recurrence:
%     T_0(x) = 1
%     T_1(x) = x
%     T_n(x) = 2*x*T_{n-1}(x) - T_{n-2}(x)

if n == 0
    T = 1;
elseif n == 1
    T = x;
else
    T0 = 1;
    T1 = x;
    for k = 2:n
        T  = 2*x*T1 - T0;
        T0 = T1;
        T1 = T;
    end
end
