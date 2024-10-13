%% Data from experimentation
Z_exp = dev6860.imps.sample{1, 2}.absz;
phase_exp = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

% Global variable used to plot the error
global error_values;
error_values = [];

% Initial values of the components
Rs = 170.2;
Rc = 6e6;
Cd = 1e-8;


x0 = [Rs Rc Cd];                     

%% Using of 'fmincon' for the optimisation of a non-linear multidimensional 
% problem with constraints 

lb=[0 0 0] ;                 % Definition of the lower bounds
ub=[inf inf inf] ;                 % upper bounds
A=[]; B=[];  Aeq=[]; Beq=[];         % No linear constraint

          
% Options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'OutputFcn', @outfun, ...
    'Algorithm', 'interior-point');

[x_opt,rmse_opt,exitflag,output] = fmincon(@(x) objective_function(x,f, ...
    Z_exp),x0,A,B,Aeq,Beq,lb,ub,@(x) constraints(x), options)

% Optimized values of the components
Rs_opt = x_opt(1);
Rc_opt = x_opt(2);
Cd_opt = x_opt(3);

fprintf('Rs optimal : %.2f\n', Rs_opt);
fprintf('Rc optimal : %.2f\n', Rc_opt);
fprintf('Cd optimal : %.2e\n', Cd_opt);
fprintf('RMSE minimal : %.4f\n', rmse_opt);

%% Plots
close all;

figure('Position', [0, 50, 600, 400]); % Impedance
Z_t = compute_abs_Zt(f, Rs, Rc, Cd);
Zt_opt = compute_abs_Zt(f, Rs_opt, Rc_opt, Cd_opt);
loglog(f,Z_exp,f,Zt_opt);
title('abs');
legend('Z-exp','Zt-opt');

figure('Position', [700, 50, 600, 400]) % Phase
phase = compute_phase(f, Rs, Rc, Cd);
phase_opt = compute_phase(f, Rs_opt, Rc_opt, Cd_opt);
semilogx(f,phase_exp,f,phase_opt);
title('phase');
legend('phase-exp','phase-opt');

figure('Position', [400, 450, 500, 300]) % Error
plot(error_values);
title('error');
xlabel('iterations');
 
%% Output function to display the error 
function stop = outfun(x, optimValues, state)
    global error_values;
    stop = false;

    if strcmp(state, 'iter')
        error_values = [error_values; optimValues.fval];
    end
end

