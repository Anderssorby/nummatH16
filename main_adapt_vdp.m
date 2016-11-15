
 %----------------------------------------------------------------------
 % Adaptive stepsize solver
 %----------------------------------------------------------------------

 
 clc
clear all 
close all

fun = @(t, y) vanDerPol(t,y);
jac = @(t,y) vdp_jac(t,y);
y0 = [2;0];
t = [0,1];
Tol = 1E-5;
 h0 = 0.1;
 
 
    % Adaptive stepsize solver
   [t,y,iflag,nfun,njac] = RKs(fun, jac, t, y0, Tol, h0);
   
   [tanal, yanal] = ode15s(fun,t,y0);
    
   err = norm(yanal(end,:)'- y(:,end));