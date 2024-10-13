function [c, ceq] = constraints(x)
    c = -x;   % x > 0
    ceq = []; % No nonlinear equality constraints
end