%% Data from experimentation
Z_exp = dev6860.imps.sample{1, 2}.absz;
phase_exp = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

% Initial values of the components
Rs = 175;
Rc = 6e6;
Cd = 5e-8;
x0 = [Rs, Rc, Cd];

objective_function(x0, f, Z_exp)

close all;
figure;
Z_t = compute_abs_Zt(f, Rs, Rc, Cd);

loglog(f,Z_exp,f,Z_t);
title('abs');
legend('Z-exp','Zt-opt');