%% Raw data from experimentation 
filename = "Data.xlsx";
data = readtable(filename);
w = 2*pi*data.freq;
mag_raw = data.mag_imp;
phase_raw = data.phase_imp;

%% export data from sweep_completo
% mag_raw = dev6860.imps.sample{1, 2}.absz;
% phase_raw = dev6860.imps.sample{1, 2}.phasez;
% f = dev6860.imps.sample{1, 2}.frequency;
% w = 2*pi*f;

[rc_fit,rs_fit,cdl_fit,alpha,tau_fit,relative_error_mag,relative_error_phase,f,...
    mag_raw,phase_raw,mag_fit,phase_fit] = algo_to_fit(w,mag_raw,phase_raw,true,0.05);

%% Different functions to adjust the graphs (optimisation)

%% Objective function used in lsqnonlin
function residuals = objective_function(x,w,y)
    a1 = x(1); % Rc (in LF domain)
    d1 = x(2); % pole : 1/Tau = 1/(Rc*Cd)
    a2 = x(3); % Rs (in HF domain)
    alpha = x(4); % Coefficient that modifies the slope
    
    int = 1 + (d1 * w * 1i).^alpha;
    residuals = y - abs(a1 ./ int + a2); % returns a vector of residuals

end 


%% A function to compute the mag and the phase using the parameters 
function [mag, phase] = compute_mag_phase(w,x)
    mag0 = zeros(length(w),1);
    a1 = x(1);
    d1 = x(2);
    a2 = x(3);
    alpha = x(4);
    
    for w_i = 1:length(w)
        int = 1 + (d1*w_i*1i).^alpha;
        mag0(w_i) = a1./int + a2;
    end 
    phase = angle(mag0); % returns the phase angle in the interval [-π,π]
    mag = abs(mag0);
end

%% Score function : return the error between raw and fit data (mag & phase)
function [error_mag_fit, error_phase_fit] = score(mag_fit,phase_fit,mag_raw,phase_raw)
    if length(mag_fit)~=length(mag_raw) || length(phase_fit)~=length(phase_raw)
        error("Can't compare, spectra of different sizes");
    end 
    error_mag_fit = abs(mean((mag_fit - mag_raw) ./ mag_raw));
    error_phase_fit = abs(mean((phase_fit - phase_raw)./ phase_raw));
end

%% Cumpute randomly the parameters 
function [lower_bound0, upper_bound0, total_init_param] = generate_random_init_param(init_param, ...
    lower_bound, upper_bound, num_init_param, absolute_lower_bound, ...
    absolute_upper_bound, var, var_bound)

    % random numbers between 0 and  var_bound (=0.05) to slightly change
    % the parameters bound
    e_pos = var_bound * rand; 
    e_neg = -var_bound * rand;
    
    % Lower bounds
    % Lower bound for the mag in LF domain 
    if lower_bound(1) * (1 + e_neg) > absolute_lower_bound(1)
        lower_bound0(1) = lower_bound(1) * (1 + e_neg);
    else
        lower_bound0(1) = lower_bound(1);
    end
    % Lower bound for the pole 
    lower_bound0(2) = 0;
    % Lower bound for the mag in HF domain 
    if lower_bound(3) * (1 + e_neg) > absolute_lower_bound(3)
        lower_bound0(3) = lower_bound(3) * (1 + e_neg);
    else
        lower_bound0(3) = lower_bound(3);
    end
    % Lower bound for alpha 
    if lower_bound(4) * (1 + e_neg) > absolute_lower_bound(4)
        lower_bound0(4) = lower_bound(4) * (1 + e_neg);
    else
        lower_bound0(4) = lower_bound(4);
    end

    % Upper bounds 
    upper_bound0(1) = Inf; % Upper bound for the mag in DC
    % Upper bound for the pole 
    if upper_bound(2) * (1 + e_pos) > absolute_upper_bound(2)
        upper_bound0(2) = upper_bound(2) * (1 + e_pos);
    else
        upper_bound0(2) = upper_bound(2);
    end
    % Upper bound for the mag in HF domain 
    if upper_bound(3) * (1 + e_pos) > absolute_upper_bound(3)
        upper_bound0(3) = upper_bound(3) * (1 + e_pos);
    else
        upper_bound0(3) = upper_bound(3);
    end
    % Upper bound for alpha
    upper_bound0(4) = 1; 

    % Computing the parameters using the bounds 
    % Pre-allocate p_i_t as a matrix
    p_i_t = zeros(length(init_param), num_init_param);
    
    for ip = 1:length(init_param)
        p = init_param(ip);
        
        % List containing random values between -var and var of size num_init_param
        e = (-var) + (2 * var) * rand(1, num_init_param);
        
        % Temporary list to store the adjusted p values
        temp_list = zeros(1, num_init_param);
        
        for i = 1:num_init_param
            p_adjusted = p * (1 + e(i));
            if upper_bound0(ip) <= p_adjusted && p_adjusted >= lower_bound0(ip)
                temp_list(i) = p_adjusted;
            else
                temp_list(i) = p;
            end
        end
        
        % Add the results to p_i_t
        p_i_t(ip, :) = temp_list;
    end
    
    total_init_param = p_i_t.' % transpose of p_i_t, contains several lists of values of
    % parameters
    
end 

%% Optimisation algorithm (main function) 
function [rc_fit,rs_fit,cdl_fit,alpha,pole_fit,relative_error_mag, relative_error_phase,...
    f,mag_raw,phase_raw,mag_fit,phase_fit] = algo_to_fit(w,mag_raw,phase_raw,plot, error_max)

    % Extract the necessary values
    mag_dc = mag_raw(1); % 1st élément of mag_raw
    pole_init = w(1) / (4 * pi); % Pole Initialisation
    mag_inf = mag_raw(end); % Last element of mag_raw
    alpha_init = 0.7; % Initial value of alpha
    
    % Initial parameters for the fit
    x0 = [mag_dc, pole_init, mag_inf, alpha_init];
    
    % Determine the lower bounds
    lower_bound = [mag_dc * 0.9,...  % Lower bound for the magnitude in DC
                      0,    ...         % Lower bound for the pôle
                      mag_inf * 0.99,... % Lower bound for the magnitude in inf
                      0.65];         % Lower bound for alpha
    
    % Determine the upper bounds
    upper_bound = [Inf,              ...      % Upper bound for the magnitude in DC
                      (w(1) / (2 * pi)) * 10,... % Upper bound for the pôle
                      mag_inf * 1.01,       ...   % Upper bound for the magnitude in inf
                      1];                     % Upper bound for alpha
    
    % Assignation of the absolute bounds
    lower_bound_abs = lower_bound;
    upper_bound_abs = upper_bound;


    
    % Options for lsqnonlin
    options = optimoptions('lsqnonlin', ...
        'Algorithm','trust-region-reflective', ...
        'TolFun', 1e-8, ...   % Tolerance on the function 
        'TolX', 1e-8, ...     % Tolerance on the gradient 
        'Display', 'iter', ...% Do not display the iterations 
        'MaxFunctionEvaluations', 2000); % Increase this from 400 to 2000     
    
    % Execute the optimisation
    [res_lsq, ~, ~, exitflag, ~] = lsqnonlin(@(x)objective_function(x,w,mag_raw), ...
        x0, lower_bound, upper_bound, options);

    [mag_fit, phase_fit] = compute_mag_phase(w, res_lsq);

    f = w/(2*pi);

    %% GENETIC ALGORITHM ( modify the parameters and redo the optimisation
    % until we reach the max number of iteration or until we reach the 
    % expected relative error )

    % Extract the values of res_lsq
    mag_dc = res_lsq(1);
    pole_fit = res_lsq(2);
    mag_inf = res_lsq(3);
    alpha_init = res_lsq(4);
    
    % Compute the adjusted values 
    rc_fit = mag_dc - mag_inf;
    rs_fit = mag_inf;
    cdl_fit = 1 / (rc_fit * pole_fit);
    alpha = alpha_init;
    
    
    % Compute the relative errors
    [relative_error_mag, relative_error_phase] = score(mag_fit, phase_fit, mag_raw, phase_raw)
    
    % Display the initial relative error
    fprintf('Initial relative error: %.2f %%\n', abs(relative_error_mag) * 100);
    
    % Initialise the parameters for the iteration
    iter_max = 5;
    init_param = x0;
    i = 0;
    
    while i < iter_max && relative_error_mag > 0.01 % relative error > 1%
        [lower_bound_i, upper_bound_i, init_param_iteration] = generate_random_init_param(init_param, ...
            lower_bound, upper_bound, 10, lower_bound_abs, upper_bound_abs, error_max, error_max);

        lower_bound_child = {};
        upper_bound_child = {};
        error_mag_child = {};

        % Loop on init_param_iteration
        for i = 1:length(init_param_iteration)
            init_param_i = init_param_iteration(i); 

            [res_lsq, ~, ~, ~, ~] = lsqnonlin(@(x)objective_function(x,w,mag_raw), x0, ...
            lower_bound_i, upper_bound_i, options);

            [mag_fit, phase_fit] = compute_mag_phase(w, res_lsq);
            
            [error_mag_fit, error_phase_fit] = score(mag_fit, phase_fit, mag_raw, phase_raw);

            % Ajouter les erreurs et les bornes aux listes
            error_mag_child{end+1} = error_mag_fit;
            lower_bound_child{end+1} = lower_bound_i;
            upper_bound_child{end+1} = upper_bound_i;

            % Mise à jour des paramètres si l'erreur est améliorée
            if relative_error_mag >= error_mag_fit
                mag_dc = res_lsq(1);
                pole_fit = res_lsq(2);
                mag_inf = res_lsq(3);
                alpha = res_lsq(4);
                
                rc_fit = mag_dc - mag_inf;
                rs_fit = mag_inf;
                cdl_fit = 1 / (rc_fit * pole_fit);
                
                % Mettre à jour les paramètres
                
                init_param = init_param_i;
                
                % Mettre à jour les erreurs relatives
                relative_error_mag = error_mag_fit;
                relative_error_phase = error_phase_fit;
            end
        end
        
        % Sort the errors in descending order and select the top 10
        error_mag_child_numeric = cell2mat(error_mag_child); % convert to a numeric variable to use sort()
        [~, index] = sort(error_mag_child_numeric, 'descend');
        best_index = index(1:10);
        
        % Initialize arrays to store filtered bounds
        best_lower_bound = cell(1, length(best_index));
        best_upper_bound = cell(1, length(best_index));
        
        % Filter bounds based on top indices
        for i = 1:length(best_index)
            idx = best_index(i);
            best_lower_bound{i} = lower_bound_child{idx};
            best_upper_bound{i} = upper_bound_child{idx};
        end
        
        % Calculate the mean for each parameter in the filtered bounds
        best_upper_bound_matrix = cell2mat(best_upper_bound');
        best_lower_bound_matrix = cell2mat(best_lower_bound');
        lower_bound = mean(best_lower_bound_matrix, 1);
        upper_bound = mean(best_upper_bound_matrix, 1);

        % Increment i
        i = i + 1;

    end
    
    % END OF THE GENETIC ALGORITHM

    %% Plots 
    close;
    if plot
        % Display the adjust values
        fprintf('\n-----------VALUES FIT-----------------------\n');
        fprintf('Rs_fit: %.1E\n', rs_fit);
        fprintf('Rc_fit: %.1E\n', rc_fit);
        fprintf('Cd_fit: %.1E\n', cdl_fit);
        fprintf('alpha: %.1E\n', alpha_init);
        fprintf('Tau: %.1E\n', pole_fit);
        
        % Display the errors
        fprintf('-----------ERRORS FIT-----------------------\n');
        fprintf('Error mag: %.1E %%\n', abs(relative_error_mag) * 100);
        fprintf('Error phase: %.1E %%\n', abs(relative_error_phase) * 100);

        w_plot = w;
        f_plot = w_plot / (2 * pi);
        
       
        figure;
        subplot(2, 1, 1); 
        loglog(f, mag_raw, '.', 'DisplayName', 'Magnitud raw'); 
        hold on;
        loglog(f_plot, mag_fit, 'DisplayName', 'Magnitud fit'); 
        ylabel('Magnitud (Nat)');
        legend('show');
        
        subplot(2, 1, 2); 
        semilogx(f, phase_raw, '.', 'DisplayName', 'Phase raw'); 
        hold on;
        semilogx(f_plot, phase_fit, 'DisplayName', 'Phase fit'); 
        ylabel('Phase (deg)');
        xlabel('Freq (Hz)');
        legend('show');
        
        
        fprintf('Error relativo final: %.1f %%\n', abs(relative_error_mag) * 100);
        fprintf('\n');
    end

end


