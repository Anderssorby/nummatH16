% Solving Robertson problem
clear all 
close all
format long 
hav
fun = @(t, y) robertson(t,y);
jac = @(t,y) robertson_jac(t,y);
y0 = [1;0;0];
t = [0,40];

% Solving for several stepsizes 
h = [0.01, 0.001];
h_len = length(h);

 
max_err = zeros(h_len,1);
max_errmeth = zeros(h_len,1);
 % Constant stepsize solver of third order
for j = 1: h_len 
    % Constant stepsize solver
    [Ya, le, tn] = ODE23_solver(y0, t, h(j), jac, fun);

    [tanal, yanal] = ode15s(fun,t,y0);
    
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
title('Loglogerror advancing method for the Robinson problem');
hold on; 
subplot(2,1,2);loglog(h,max_errmeth,'r');
hold on; axis([min(h)/10 max(h)*10 min(max_errmeth)/10 max(max_errmeth)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror error method for the Robinson problem')
