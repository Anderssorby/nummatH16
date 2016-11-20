system()

syms z gamma yn

% lambda*h = z

Y1 = yn;

Y2 = 1/(1-z*gamma)*(yn + z*a(2,1)*Y1);

Y3 = 1/(1-z*gamma)*(yn + z*a(3,1)*Y1 + z*a(3,2)*Y2);

Y4 = 1/(1-z*gamma)*(yn + z*a(4,1)*Y1 + z*a(4,2)*Y2 + z*a(4,3)*Y3);

simplify(Y3/yn)

simplify(Y4/yn)

