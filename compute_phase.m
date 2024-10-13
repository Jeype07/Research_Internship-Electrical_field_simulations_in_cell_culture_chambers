function phase = cumpute_phase(f, Rs, Rc, Cd)
  % Compute the abs given some parameters
  w = 2 * pi * f;
  tau = Rc * Cd;
  phase = atan2(-Rc*tau*w./(1+tau^2*w.^2), (Rs + Rc +tau^2*w.^2*Rs)./(1+tau^2*w.^2));
end
