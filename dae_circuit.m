function [ yd ] = dae_circuit(t, y, V, R, A,C)



yd = (V(t) - y)*(1+A)/(R*C);



end