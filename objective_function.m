function rmse_value = objective_function(x, f, exp_data)
    Rs = x(1);
    Rc = x(2);
    Cd = x(3);
    format long e;
    tau = Rc * Cd;
    w = 2 * pi * f;
    % Optimisation using the impedance's formula 
    Z_t = sqrt(((Rs + Rc + tau^2*w.^2*Rs)./(1+tau^2*w.^2)).^2 + ...
        ((Rc*tau*w)./(1+tau^2*w.^2)).^2);
    rmse_value = sqrt(mean((log(Z_t) - log(exp_data)).^2));

    % Optimisation using the phase's formula
%     phase = atan2(-Rc*tau*w./(1+tau^2*w.^2), ...
%         (Rs + Rc +tau^2*w.^2*Rs)./(1+tau^2*w.^2));
%     rmse_value = sqrt(mean((phase - exp_data).^2));
end

