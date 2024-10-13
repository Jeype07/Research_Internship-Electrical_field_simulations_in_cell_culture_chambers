z_real = dev6860.imps.sample{1, 2}.absz;
z_img = dev6860.imps.sample{1, 2}.phasez;
f = dev6860.imps.sample{1, 2}.frequency;

figure(1);
loglog(f,z_real);

figure(2);
semilogx(f,z_img);