
function t = atan2d_0_360(y,x)
   if x == 0
       if y == 0
           t = 0;
       elseif y > 0
           t = 90;
       else
           t = 270;
       end
    elseif x > 0
       if y >= 0
           t = atand(y/x);
       else
           t = atand(y/x) + 360;
       end
    elseif x < 0
       if y == 0
           t = 180;
       else
           t = atand(y/x) + 180;
       end
    end
end