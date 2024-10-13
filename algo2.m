% Data from experimentation
Z_exp = dev6860.imps.sample{1, 2}.absz;
phase = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;


% Initial values of the components
Rs = 175;
Rc = 6e6;
Cd = 5e-8;
x0 = [Rs, Rc, Cd];

% Fonction objectif à minimiser
objective = @(x) rmse_loss_log(Zt(f, x(1), x(2), x(3)), Z_exp);

% Limits for the values of the parameters
lb = [1e2, 1e6, 1e-9];
ub = [1e3, 1e7, 1e-8];

% Optimisation
%x_opt = lsqnonlin(objective, x0, lb, ub);
options = optimoptions(@particleswarm, 'Display', 'iter', 'MaxIterations', 100);
[x_opt, fval] = particleswarm(objective, 3, lb, ub, options)

% % Extract the optimized parameters
Rs_opt = x_opt(1)
Rc_opt = x_opt(2)
Cd_opt = x_opt(3)


% Display the curves
close all;
figure;
Z_t = Zt(f, Rs, Rc, Cd);
Zt_opt = Zt(f, Rs_opt, Rc_opt, Cd_opt);
loglog(f,Z_exp,f,Zt_opt);
title('abs');
legend('Z-exp','Zt-opt');

% figure;
% phase = compute_phase(f, Rs, Rc, Cd);
% phase_opt = compute_phase(f, Rs_opt, Rc_opt, Cd_opt);
% semilogx(f,phase_exp,f,phase_opt);
% title('phase');
% legend('phase-exp','phase-opt');


% Function of the electrical model of Z
function Z = Zt(f, Rs, Rc, Cd)
  w = 2 * pi * f;
  tau = Rc * Cd;
  Z= sqrt(((Rs + Rc +tau^2*w.^2*Rs)./(1+tau^2*w.^2)).^2 + ((Rc*tau*w)./(1+tau^2*w.^2)).^2);
end
