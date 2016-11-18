% Solving Van der Pol problem
clc
clear all 
close all

my = 50;
fun = @(t, y) vanDerPol(t,y, my);
jac = @(t,y) vdp_jac(t,y, my);
y0 = [1;2];
t = [0,my];
Tol = 1E-5;

% Solving for several stepsizes 
h = [0.1, 0.01, 0.001, 0.0001, 0.00001];
h_len = length(h);

 
max_err = zeros(h_len,1);
max_errmeth = zeros(h_len,1);
 % Constant stepsize solver of third order
 
     options = odeset('RelTol',1e-8,'AbsTol',1e-10);
    [tanal, Yanal] = ode15s(fun,t,y0, options);
for j = 1: h_len 
    % Constant stepsize solver
    [Ya, le, Y3, tn] = ODE23_solver(y0, t, h(j), jac, fun);
    [Y3, tn] = ODE2_solver(y0, t, h(j), jac, fun);

    
    
    max_err(j) = max(norm(Yanal(end,:)'-Ya(:,end)));
    max_errmeth(j) = max(norm(Yanal(end,:)'-Y3(:,end)));
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
title('Loglogerror advancing method for VDP problem');

%Add line to show order of method

line = zeros(2,1);
line(1) = h(1)^3/10000;
line(2) = h(end)^3/10000;
hvec = [h(1);h(end)];


loglog(hvec, line);

%Third order plot


hold on; 
subplot(2,1,2);loglog(h,max_errmeth,'r');
hold on; axis([h(h_len)/10 h(1)*10 min(max_errmeth)/10 max(max_errmeth)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror error method for linear test problem')

%Add line to show order of method


line = zeros(2,1);
line(1) = h(1)^2/1000;
line(2) = h(end)^2/1000;
hvec = [h(1);h(end)];


loglog(hvec, line);


 %----------------------------------------------------------------------
 % Adaptive stepsize solver
 %----------------------------------------------------------------------

 h0 = 0.1;
 
    % Adaptive stepsize solver
  %  [t,y,iflag,nfun,njac] = RKs(fun, jac, t, y0, Tol, h0);
    
 %   err = norm(yanal(end,:)'- y(:,end))