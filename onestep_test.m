%x7
clc;
clear all
format long 

% Correctness linear test function

     Tolit = 1e-6;
     jac = [-2 1; 1 -2];
     h = 0.1;
     tn = 0;
     yn = [1;2];
     
     f = @(t, y) lintest(t,y);
     
     [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit);
     fprintf('ynext = %f\n', ynext);
     fprintf('tnext = %f\n', tnext);
     yAnal = analLinTest(tnext)
     
     
% Test van der pool equation

     
     h = 0.1;
     tn = 0;
     yn = [2;0];
     my = 50;
     
     jac = [0,1;-my*yn(1)*yn(2)-1, my*(1-yn(1)^2)];
     
     v = @(t,y) vanDerPol(t,y, my);
     
     [tnext, ynext, le, iflag] = onestep(v, jac, tn, yn, h, Tolit);
     fprintf('ynext = %f\n', ynext);
     fprintf('tnext = %f\n', tnext);