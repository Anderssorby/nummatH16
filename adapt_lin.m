% Solving linear test problem with adaptive method

clear all
close all

fun = @(t, y) lintest(t,y);
jac = @(t, y) [-2 1; 1 -2];
y0 = [1;2];
t = [0,1];
tend = t(2);
Tol = 1e-2;

 h0 = 0.1;%t(2)-t(1);
 
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs(fun, jac, t, y0, Tol, h0);
    
    %Petesolution(Does not werk with our code)
    %[t,y,iflag,nfun,njac] = RKsFromPete(fun, jac, 0,1, y0, Tol, h0);
    
    
    % Analytic solution
    Yanal = analLinTest(t(2));
    
    % Work
    werk = nfun + 2*njac
    
    % End time stuff
    endTime = t(end)
    diffEndTime = endTime - tend
    
    
    % Error
    err = norm(Yanal- y(:,end))