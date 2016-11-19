% Solving Van der Pol problem
clc
clear all 
close all

my = 50;
fun = @(t, y) vanDerPol(t,y, my);
jac = @(t,y) vdp_jac(t,y, my);
y0 = [2;0];
t = [0,1];

% Solving for several stepsizes 
h = [0.1 0.01 0.001 0.0001];
h_len = length(h);

 
max_err = zeros(h_len,1);
max_errmeth = zeros(h_len,1);
 % Constant stepsize solver of third order
    options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [tanal, yanal] = ode15s(fun,t,y0, options);
for j = 1: h_len 
    % Constant stepsize solver
    [Ya, le, tn] = ODE23_solver(y0, t, h(j), jac, fun);

    Ya(:,end)
    max_err(j) = norm(yanal(end,:)'-Ya(:,end));
    max_errmeth(j) = max(le(:,end));
end


% Checking that the order conditions are right 
ord = polyfit(log(h),log(max_err'), 1);
pnum = ord(1)
ord = polyfit(log(h),log(max_errmeth'), 1);
pnum = ord(1)

% Creating convergence plot
figure; 
subplot(2,1,1);loglog(h,max_err);
hold on; axis([h(h_len)/10 h(1)*10 min(max_err)/10 max(max_err)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror advancing method for the van der pol problem');
hold on; 
subplot(2,1,2);loglog(h,max_errmeth,'r');
hold on; axis([min(h)/10 max(h)*10 min(max_errmeth)/10 max(max_errmeth)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror error method for the van der pol problem')



 %----------------------------------------------------------------------
 % Adaptive stepsize solver
 %----------------------------------------------------------------------

 h0 = 0.1;
 
    % Adaptive stepsize solver
  %  [t,y,iflag,nfun,njac] = RKs(fun, jac, t, y0, Tol, h0);
    
 %   err = norm(yanal(end,:)'- y(:,end))