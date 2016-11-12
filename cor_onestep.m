%x7
clc;
clear all
format long 

% Correctness linear test function

     Tolit = 1e-6;
     jac = @(t,y) [-2 1; 1 -2];
     h = 1;
     tn = 0;
     yn = [1;2];
     
     f = @(t, y) lintest(t,y);
     
     [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit);
     fprintf('ynext = %f\n', ynext);
     fprintf('tnext = %f\n', tnext);
     yAnal = analLinTest(tnext)
     
     
fprintf('Test van der pool equation')

     
     h = 0.1;
     tn = 0;
     yn = [2;0];
     my = 50;
     
     jac = @(t,y) [0,1;-my*yn(1)*yn(2)-1, my*(1-yn(1)^2)];
     
     v = @(t,y) vanDerPol(t,y, my);
     
     [tnext, ynext, le, iflag] = onestep(v, jac, tn, yn, h, Tolit);
     fprintf('ynext = %f\n', ynext);
     fprintf('tnext = %f\n', tnext);
     
     
fprintf('Test Robertson reaction')

     
     h = 0.1;
     tn = 0;
     yn = [1;0;0];
     
     jac = @(t,y) [
        (-0.04 + 10^4*y(2)*y(3)) (-0.04*y(1) + 10^4*y(3)) (-0.04*y(1) + 10^4*y(2));
        (0.04 - 10^4*y(2)*y(3) - 3*10^7*y(2)^2)  (0.04*y(1) - 10^4*y(3) - 6*10^7*y(2)) (0.04*y(1) - 10^4*y(2) - 3*10^7*y(2)^2);
        0 6*10^7*y(2) 0;
     ];
     
     f = robertson;
     
     [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit);
     fprintf('ynext = %f\n', ynext);
     fprintf('tnext = %f\n', tnext);
     
  
