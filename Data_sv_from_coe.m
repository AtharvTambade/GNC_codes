deg = pi/180;
mu = 398600;
%...Data declaration for Example 4.5 (angles in degrees):
h = 52884;
e = 0.0074922;
RA = 127.4561*deg; %Right ascension of the node (radians)
i = 64.5538*deg; %Inclination (radians)
w = 239.1985*deg; %Argument of perigee (radians)
TA = 2.1103; %True anomaly (radians)
%...
coe = [h, e, RA, i, w, TA];
%...Algorithm 4.5 (requires angular elements be in radians):
[r, v] = sv_from_coe(coe, mu);
%...Echo the input data and output the results to the command window:
fprintf('Gravitational parameter (km^3/s^2) = %g\n', mu)
fprintf('Angular momentum (km^2/s) = %g\n', h)
fprintf('Eccentricity = %g\n', e)
fprintf('Right ascension (deg) = %g\n', RA)
fprintf('Argument of perigee (deg) = %g\n', w)
fprintf('True anomaly (deg) = %g\n', TA)
fprintf('State vector:\n')
fprintf('r (km) = [%g %g %g]\n', r(1), r(2), r(3))
fprintf('v (km/s) = [%g %g %g]\n', v(1), v(2), v(3))



                                                    