%% Data from experimentation
Z_exp = dev6860.imps.sample{1, 2}.absz;
phase_exp = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

% Global variable used to plot the error
global error_values;
error_values = [];

% Initialisation of variables and values of the components
List_Rs = linspace(150,200,10);
List_Rc = linspace(1e6,1e7,10);
List_Cd = linspace(1e-8,1e-6,10);
format long e;

rmse_value = objective_function([List_Rs(1) List_Rc(1) List_Cd(1)],f,phase_exp);
Comp_values = [List_Rs(1) List_Rc(1) List_Cd(1)];


%% Using of 'fmincon' for the optimisation of a non-linear multidimensional problem with constraints 
for i = 1:length(List_Rs)
    i
    for j = 1:length(List_Rc)
        for k = 1:length(List_Cd)
            
            Rs = List_Rs(i);
            Rc = List_Rc(j);
            Cd = List_Cd(k);

            x0 = [Rs Rc Cd];                     
            
            
            
            lb=[List_Rs(1) List_Rc(1) List_Cd(1)];                 % Definition of the lower bounds
            ub=[List_Rs(length(List_Rs)) List_Rc(length(List_Rc)) List_Cd(length(List_Cd))];            % upper bounds
            A=[]; B=[];  Aeq=[]; Beq=[];         % No linear constraint
            
            % Validate the initial point
            if any(x0 < lb) || any(x0 > ub)
                error('Initial point x0 is not within the bounds lb and ub.');
            end


            % Options for fmincon
            options = optimoptions('fmincon', 'Display', 'none');
            
            [x_opt,rmse_opt,exitflag,output] = fmincon(@(x) objective_function(x,f,Z_exp),x0,A,B,Aeq,Beq,lb,ub,@(x) constraints(x),options);
            
            if rmse_opt < rmse_value
                Comp_values = x0;
                rmse_value = rmse_opt

            end


        end
    end
end
% Optimized values of the components
Rs_opt = Comp_values(1);
Rc_opt = Comp_values(2);
Cd_opt = Comp_values(3);

fprintf('Rs optimal : %.2f\n', Rs_opt);
fprintf('Rc optimal : %.2f\n', Rc_opt);
fprintf('Cd optimal : %.2e\n', Cd_opt);
fprintf('RMSE minimal : %.4f\n', rmse_value);

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

 
% %% Output function to display the error 
% function stop = outfun(x, optimValues, state)
%     global error_values;
%     stop = false;
% 
%     if strcmp(state, 'iter')
%         error_values = [error_values; optimValues.fval];
%     end
% end

