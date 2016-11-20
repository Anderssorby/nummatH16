function [ jac ] = dae_circuit_jac(t, y, R, A,C)


jac = -(1+A)/(R*C);


end