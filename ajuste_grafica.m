abs = dev6860.imps.sample{1, 2}.absz;
phase = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

Rs = 175;
Rc = [20e6, 23e6, 26e6, 30e6];
Cd = 50e-9;


w=2*pi*f;

abs_fitted = zeros(4,300);
phase_fitted = zeros(4,300);

for i=1:4
    tau = Cd*Rc(i);
    abs_fitted(i,:) = sqrt(((Rs + Rc(i) +tau^2*w.^2*Rs)./(1+tau^2*w.^2)).^2 + ((Rc(i)*tau*w)./(1+tau^2*w.^2)).^2);
    phase_fitted(i,:) = atan2(-Rc(i)*tau*w./(1+tau^2*w.^2), (Rs + Rc(i) +tau^2*w.^2*Rs)./(1+tau^2*w.^2));
end

figure(1);
loglog(f,abs);
hold on;
loglog(f,abs_fitted(:,:));
legend(["Measured","20Meg","23Meg","26Meg","30Meg"]);

figure(2);
semilogx(f,phase);
hold on;
semilogx(f,phase_fitted(:,:));
legend(["Measured","20Meg","23Meg","26Meg","30Meg"]);