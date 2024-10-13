abs = dev6860.imps.sample{1, 2}.absz;
phase = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

Rs = 175;
Rc = 23e6;
Cd = 50e-9;

tau = Cd*Rc;
w=2*pi*f;

abs_fitted = sqrt(((Rs + Rc +tau^2*w.^2*Rs)./(1+tau^2*w.^2)).^2 + ((Rc*tau*w)./(1+tau^2*w.^2)).^2);
phase_fitted = atan2(-Rc*tau*w./(1+tau^2*w.^2), (Rs + Rc +tau^2*w.^2*Rs)./(1+tau^2*w.^2));


figure(1);
loglog(f,abs);
hold on;
loglog(f,abs_fitted);

figure(2);
semilogx(f,phase);
hold on;
semilogx(f,phase_fitted);