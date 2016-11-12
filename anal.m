syms x(t) y(t)

A = [-2 1;1 -2];

B = [t; t+3];

Y = [x;y];

eqn = diff(Y) == A*Y+B;

[xsol(t) ysol(t)]  = dsolve(eqn);


exp(-t)*(C2 + (exp(t)*(2*t + 1))/2) - exp(-3*t)*(C1 + exp(3*t)/2);
exp(-t)*(C2 + (exp(t)*(2*t + 1))/2) + exp(-3*t)*(C1 + exp(3*t)/2);
