%% Main Robertson reaction
% ------------------------------------------------------------------------
% Solving the Robertson reaction with adaptive step size

clear all
close all


%% Variable definitions
fun = @(t, y) robertson(t,y);
jac = @(t, y) robertson_jac(t,y);
y0 = [1;0;0]; % Initial conditions. 
tint = [0,40]; % Time interval
Tol = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
len_Tol = length(Tol); 
h0 = 5; % Initial first step 


%% Solver
for i = 1:len_Tol
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs(fun, jac, tint, y0, Tol(i), h0);
   
    
    % Stiff built-in Matlab solver ODE15s.
     options = odeset('RelTol',1e-8,'AbsTol',1e-10);
     [tanal, yanal] = ode15s(fun,tint,y0, options);
    
    % Work
    werk(i) = nfun + 3*njac;
    
    % Error 
    err(i) = norm(yanal(end,:)'- y(:,end));
    
end
      
    
%% Postprocess
    
    figure(); loglog(Tol, err); grid on
    hold on 
    xlabel('Tolerance'); ylabel('Error');
    title('Error Robertson reaction');    
    
    figure(); loglog(Tol , werk); grid on 
    hold on
    xlabel('Tolerance'); ylabel('Work');
    title('Work Robertson reaction');
    
    figure(); plot(t,y(1,:), '-o'); 
    hold on 
    xlabel('t');ylabel('y1');
    title('Robertson reaction'); 
    
    figure(); plot(t,y(2,:), '-o',tanal,yanal(:,2),'r');
    hold on
    xlabel('t'); ylabel('y2');
    title('Robertson reaction');
    legend('ODE23','ODE15s');
    
    
    figure();plot(t, y(1,:), 'b');
    hold on 
    xlabel('t'); ylabel('y1');
    title('Robertson reaction');
    
    figure();plot(tanal, yanal(:,1), 'b');
    hold on 
    xlabel('t'); ylabel('y1');
    title('Robertson reaction');
    
    figure();plot(t, y(3,:), 'b');
    hold on 
    xlabel('t'); ylabel('y3');
    title('Robertson reaction');
    
    figure();plot(tanal, yanal(:,3), 'b');
    hold on 
    xlabel('t'); ylabel('y3');
    title('Robertson reaction');
    
    figure(); 