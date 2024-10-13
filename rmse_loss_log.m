function rmse = rmse_loss_log(f, Rs, Rc, Cd, Z_exp)
  Zt = compute_abs_Zt(f, Rs, Rc, Cd);
  
  %Calculate the RMSE
  rmse = sqrt(mean((abs(log(Zt)) - abs(log(Z_exp))).^2));
end

