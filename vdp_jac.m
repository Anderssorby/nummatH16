function jac = vdp_jac(t,y)
    my = 50;
    jac = [0,                   1;...
           -2*my*y(1)*y(2)-1,    my*(1-+
       y(1)^2)];
end