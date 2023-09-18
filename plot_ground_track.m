function plot_ground_track
% 
   hold on
   xlabel('East Longitude(degrees)')
   ylabel('Latitude (degrees)')
   axis equal
   grid on
   for i = 1:n_curves
      plot(RA{i}, Dec{i})
   end
   axis ([0 360 -90 90])
   text( ra(1), dec(1), 'o Start')
   text(ra(end), dec(end), 'o Finish')
   line([min(ra) max(ra)],[0 0], 'Color','k') %the equator
end
%plot_ground_track
% 