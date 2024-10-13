function Zt = compute_abs_Zt(f, Rs, Rc, Cd)
  % Compute the abs given some parameters
  w = 2 * pi * f;
  tau = Rc * Cd;
  Zt = sqrt(((Rs + Rc +tau^2*w.^2*Rs)./(1+tau^2*w.^2)).^2 + ((Rc*tau*w)./(1+tau^2*w.^2)).^2);
end

