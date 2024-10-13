function rmse = rmse_loss_log(Zt, Z_exp)
  % Ensure both inputs are column vectors
  Zt = Zt(:);
  Z_exp = Z_exp(:);

  % Calculate the RMSE
  rmse = sqrt(mean((abs(Zt) - abs(Z_exp)).^2));
end

