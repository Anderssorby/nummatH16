function [jac] = robertson_jac(t,y)

jac = [-0.04,               1e4*y(3),           1e4*y(2);...
        0.04,     -1e4*y(3)-6e7*y(2),          -1e4*y(2);...
           0,               6e7*y(2),                  0];

end