function print_orbital_data
% 
coe = [h e Wo incl wpo TAo];
[ro, vo] = sv_from_coe(coe, mu);
fprintf('Angular momentum = %g km^2/s\n' , h)
fprintf('Eccentricity = %g\n' , e)
fprintf('Semimajor axis = %g km\n' , a)
fprintf('Perigee radius = %g km\n' , rP)
fprintf('Apogee radius = %g km\n' , rA)
fprintf('Period = %g hours\n' , T/3600)
fprintf('Inclination = %g deg\n' , incl/deg)
fprintf('Initial true anomaly = %g deg\n' , TAo/deg)
fprintf('Time since perigee = %g hours\n' , to/3600)
fprintf('Initial RA = %g deg\n' , Wo/deg)
fprintf('RA_dot = %g deg/period\n' , Wdot/deg*T)
fprintf('Initial wp = %g deg\n' , wpo/deg)
fprintf('wp_dot = %g deg/period\n' , wpdot/deg*T)
fprintf('r0 = [%12g, %12g, %12g] (km)\n', ro(1), ro(2), ro(3))
fprintf('magnitude = %g km\n\n', norm(ro))
fprintf('v0 = [%12g, %12g, %12g] (km)\n', vo(1), vo(2), vo(3))
fprintf('magnitude = %g km\n\n', norm(vo))
end %print_orbital_data