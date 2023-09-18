function find_ra_and_dec
% 
% Propagates the orbit over the specified time interval, transforming
% the position vector into the earth-fixed frame and, from that,
% computing the right ascension and declination histories
% ––––––––––––––––––––––
%
   times = linspace(to,tf,1000);
   ra = [];
   dec = [];
   theta = 0;
   for i = 1:length(times)
      t = times(i);
      M = 2*pi/T*t;
      E = kepler_E(e, M);
      TA = 2*atan(tan(E/2)*sqrt((1+e)/(1-e)));
      r =h^2/mu/(1 + e*cos(TA))*[cos(TA) sin(TA) 0].';                                                                                                                                                                                                                                                               
      W = Wo + Wdot*t;
      wp = wpo + wpdot*t;
      R1 = [ cos(W) sin(W) 0; -sin(W) cos(W) 0; 0 0 1];
      R2 = [1 0 0; 0 cos(incl) sin(incl); 0 -sin(incl) cos(incl)];
      R3 = [ cos(wp) sin(wp) 0; -sin(wp) cos(wp) 0; 0 0 1];
      QxX = (R3*R2*R1).';
      R = QxX*r;
      theta = we*(t - to);
      Q = [ cos(theta) sin(theta) 0; -sin(theta) cos(theta) 0; 0 0 1];
      r_rel = Q*R;
      [alpha, delta] = ra_and_dec_from_r(r_rel);
      ra = [ra; alpha];
      dec = [dec; delta];
   end
end
%find_ra_and_dec
% 