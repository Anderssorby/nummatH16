%% Main Linear test problem equation
% ------------------------------------------------------------------------
% Solving the linear test equation with adaptive step size

clear all
close all


%% Variable definitions
fun = @(t, y) lintest(t,y);
jac = @(t, y) lintest_jac(t,y);
y0 = [1;2]; % Initial conditions. 
tint = [0,1]; % Time interval
Tol = [1e-2 1e-3 1e-4 1e-5 1e-6 1e-7 1e-8];
len_Tol = length(Tol); 
h0 = 0.1; % Initial first step 


%% Solver
for i = 1:len_Tol
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs(fun, jac, tint, y0, Tol(i), h0);
   
    
    %Analytic solution
    yanal = lintest_anal(t);
    
    % Work
    werk(i) = nfun + 2*njac;
    
    % Error 
    err(i) = norm(yanal(:,end)-y(:,end));
    
end
      
    
%% Postprocess
    
    figure(); loglog(Tol, err);
    hold on 
    xlabel('Tolerance'); ylabel('Error');
    title('Error linear test equation');    
    
    figure(); loglog(Tol , werk);
    hold on
    xlabel('Tolerance'); ylabel('Work');
    title('Work linear test equation');
    
    figure(); plot(t,y(1,:), '-o');
    hold on 
    xlabel('t');ylabel('y1');
    title('Linear test equation'); 
    
    figure(); plot(t,y(2,:), '-o',t,yanal(2,:),'r');
    hold on
    xlabel('t'); ylabel('y2');
    title('Linear test equation');
    legend('ODE23','Analytic');
    
    
    figure();plot(t, y(1,:), 'b');
    hold on 
    xlabel('t'); ylabel('y1');
    title('Linear test equation');