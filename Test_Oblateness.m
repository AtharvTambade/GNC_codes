r = [-2429.1 4555.1 4577];
v = [ -4.7689 -5.6113 3.0535];
hours = 72;
[r_fin, v_fin] = rv_final_from_rv_initial(r,v,hours);
disp(r_fin);
disp(v_fin);