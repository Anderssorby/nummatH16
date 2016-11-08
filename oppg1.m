system()

syms h la gamma yn

Y1 = yn

Y2 = yn/(1-h*gamma)*(1 + h*a(2,1))

Y3 = 1/(1-h*gamma)*(yn+a(3,1)*yn*h + h*a(3,2)*Y2)

Y4 = 1/(1-h*gamma)*(h*a(4,1)*yn + h*a(4,2)*Y2 + h*a(4,3)*Y3)

simplify(Y3 / yn)


simplify(Y4)

simplify(Y4/yn)
