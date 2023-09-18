function form_separate_curves
% 
% Breaks the ground track up into separate curves which start
% and terminate at right ascensions in the range [0,360 deg].
% –––––––––––––––––––––––––––
   tol = 100;
   curve_no = 1;
   n_curves = 1;
   k = 0;
   ra_prev = ra(1);
   for i = 1:length(ra)
      if abs(ra(i) - ra_prev) > tol 
         curve_no = curve_no + 1;
         n_curves = n_curves + 1;
         k = 0;
      end
      k = k + 1;
      RA{curve_no}(k) = ra(i);
      Dec{curve_no}(k) = dec(i);
      ra_prev = ra(i);
   end
end
%form_separate_curves
% 