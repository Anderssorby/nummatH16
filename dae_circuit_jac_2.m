function [ yd ] = dae_circuit_jac_2(t, y, R, A)



yd = [0 0 0 0 -1;...
      0 -1/R 1/R 0 0;...
      0 1/R -1/R 1 0;...
      0 -A -1 0 0;...
      -1 0 0 0 0];



end