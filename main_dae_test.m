% Main DAE test problem
% ------------------------------------------------------------------------
% Solving the Robertson reaction with adaptive step size

clear all
close all


%% Variable definitions
fun = @(t, y) dae_test(t,y);
jac = @(t, y) dae_test_jac(t,y);
M = [1 -1 0; -1 1 0; 0 0 0];
y0 = [0;-1;4]; % Initial conditions. 
tint = [0,1]; % Time interval
Tol = [1e-4 1e-5];
len_Tol = length(Tol); 
h0 = 0.1; % Initial first step 


%% Solver
for i = 1:len_Tol
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs_tweak(fun, jac, M, tint, y0, Tol(i), h0);
   
    
    % Stiff built-in Matlab solver ODE15s.
%     options = odeset('RelTol',1e-8,'AbsTol',1e-10);
%     [tanal, yanal] = ode15s(fun,tint,y0, options);
    
    % Work
    werk(i) = nfun + 3*njac;
    
    % Error 
    err(i) = norm(yanal(end,:)'- y(:,end));
    
end
      
figure();plot(t,y(1,:))
figure();plot(t,y(2,:))
figure();plot(t,y(3,:))