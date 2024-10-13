% Data from experimentation
Z_exp = dev6860.imps.sample{1, 2}.absz;
phase = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

% Initial values of the components
% Rs = 175;
% Rc = 23e6;
% Cd = 50e-9;
Rs = 175;
Rc = 6e6;
Cd = 5e-8;

% RMSE (objective function or function to minimise)
err = @(x, y, z) rmse_loss_log(f, x, y, z, Z_exp);
err(Rs, Rc, Cd);

% Define small increments (dRs, dRc, dCd)
dRs = 0.1;   % Increment for dRs in 10^2 Ohms
dRc = 1e3;  % Increment for dRc in MOhms
dCd = 1e-12;  % Increment for dCd in nF

% Calculate the partial derivatives using finite differences
% Partial derivative with respect to dRs
derr_dRs = @(x, y, z) (err(x + dRs, y, z) - err(x, y, z)) / dRs;

% Partial derivative with respect to dRc
derr_dRc = @(x, y, z) (err(x, y + dRc, z) - err(x, y, z)) / dRc;

% Partial derivative with respect to dCd
derr_dCd = @(x, y, z) (err(x, y, z + dCd) - err(x, y, z)) / dCd;


% Gradient descent : xk+1 = xk - h*grad(err(xk))
% Initial guess : 
comp0 = [Rs, Rc, Cd];
h = 0.01;
max_iter = 100;
List_err = zeros(1,max_iter+1);
List_err(1) = err(comp0(1), comp0(2), comp0(3));
comp_old = comp0;
for k = 1:max_iter
    comp_new = comp_old - h*[derr_dRs(comp_old(1), comp_old(2), comp_old(3)), derr_dRc(comp_old(1), comp_old(2), comp_old(3)), derr_dCd(comp_old(1), comp_old(2), comp_old(3))];
    List_err(k+1) = err(comp_new(1), comp_old(2), comp_new(3));
end 

close all;
plot(List_err, '-*', 'LineWidth', 1.5, 'MarkerSize', 5);
comp_new