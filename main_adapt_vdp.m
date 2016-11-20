%% Main Van der Pol equation
% ------------------------------------------------------------------------
% Solving the Van der Pol equation with adaptive step size

clear all
close all


%% Variable definitions
my = 50;
fun = @(t, y) vanDerPol(t,y, my);
jac = @(t, y) vdp_jac(t,y, my);
y0 = [2;0]; % Initial conditions. 
tint = [0,my]; % Time interval
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
    werk(i) = nfun + 2*njac;
    werkmatlab(i) = yanal.stats.nfevals;
    
    % Error 
    err(i) = norm(yanal.y(:,end)- y(:,end));
    
end
      
    
%% Postprocess
    
    figure(); 
    loglog(Tol, err); grid on
    hold on 
    xlabel('Tolerance'); ylabel('Error'); 
    title('Error Van der Pol oscillator, my = 50');    
    
    figure(); 
    loglog(err, werk); grid on 
    hold on
    loglog(err,werkmatlab);
    xlabel('error'); ylabel('Work');
    title('Work precision diagram Van der Pol oscillator, my = 50');
    legend('ODE23','ODE15s')
    
    figure(); 
    subplot(2,1,1)
    plot(t,y(1,:), '-o');
    hold on 
    xlabel('t');ylabel('y1');
    title('Van der Pol oscillator, my = 50'); 
    
    hold on
    subplot(2,1,2)
    plot(t,y(2,:), '-o',yanal.x,yanal.y(2,:),'r');
    hold on
    xlabel('t'); ylabel('y2');
    title('Van der Pol oscillator, my = 50');
    legend('ODE23','ODE15s');
    
    figure();
    plot(t, y(1,:), 'b');
    hold on 
    xlabel('t'); ylabel('y1');
    title('Van der Pol oscillator, my = 50');