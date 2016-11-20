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


% Stiff built-in Matlab solver ODE15s.
options = odeset('RelTol',1e-10,'AbsTol',1e-12);
[yanal] = ode15s(fun,tint,y0, options);
yanal.stats
     
%% Solver
for i = 1:len_Tol
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs(fun, jac, tint, y0, Tol(i), h0);
   
    % Work
    werk(i) = nfun + 3*njac;
    werkmatlab(i) = yanal.stats.nfevals;
    
    % Error 
    err(i) = norm(yanal.y(:,end) - y(:,end));
end
      
    
%% Postprocess
    
    figure(); 
    loglog(Tol, err); grid on
    hold on 
    xlabel('Tolerance'); ylabel('Error');
    title('Error Robertson reaction');    
    
    figure(); 
    loglog(err , werk); grid on 
    hold on
    loglog(err, werkmatlab);
    xlabel('Error'); ylabel('Work');
    title('Work precision diagram Robertson reaction');
    
    figure(); 
    subplot(3,1,1)
    plot(t,y(1,:), '-o',yanal.x,yanal.y(1,:),'r'); 
    hold on 
    xlabel('t');ylabel('y1');
    title('Robertson reaction'); 
    
    hold on 
    subplot(3,1,2)
    plot(t,y(2,:), '-o',yanal.x,yanal.y(2,:),'r');
    hold on
    xlabel('t'); ylabel('y2');
    title('Robertson reaction');
    legend('ODE23','ODE15s');
    
    hold on 
    subplot(3,1,3)
    plot(t,y(3,:), '-o',yanal.x,yanal.y(3,:),'r');
    hold on
    xlabel('t'); ylabel('y2');
    title('Robertson reaction');
    legend('ODE23','ODE15s');
 