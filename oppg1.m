rksystem()

syms h la gamma yn

Y1 = yn

Y2 = yn/(1-h*la*gamma)*(1 + h*la*a(2,1))

Y3 = 1/(1-h*la*gamma)*(yn+a(3,1)*yn*h + la*h*a(3,2)*Y2)

Y4 = 1/(1-h*la*gamma)*(la*h*a(4,1)*yn + la*h*a(4,2)*Y2 + la*h*a(4,3)*Y3)

simplify(Y3 / yn)


simplify(Y4)

simplify(Y4/yn)
