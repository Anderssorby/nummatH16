% Solving linear test problem with adaptive method

%clear all
close all


my = 50;
fun = @(t, y) vanDerPol(t,y, my);
jac = @(t, y) vdp_jac(t,y, my);
y0 = [2;0];
t = [0,my];
Tol = 1e-7;


 h0 = Tol*100;
 
    % Adaptive stepsize solver
    %[tn,y,iflag,nfun,njac] = RKs(fun, jac, t, y0, Tol, h0);
    
    %Petesolution(Does not werk with our code)
    [tn,y,iflag,nfun,njac] = RKsFromPete(fun, jac, 0,1, y0, Tol, h0);
    
    % Analytic solution
%     options = odeset('RelTol',1e-8,'AbsTol',1e-10);
%     [tanal, Yanal] = ode15s(fun,t,y0, options);
    
    % Work
    werk = nfun + 2*njac
    
    % End time stuff
    endTime = tn(end)
    diffEndTime = endTime - tanal(end)
    
    % Error 
    err = norm(Yanal(end,:)'- y(:,end))
    