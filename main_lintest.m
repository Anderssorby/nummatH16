% Solving linear test problem
clc
clear all
close all

fun = @(t, y) lintest(t,y);
jac = @(t, y) [-2 1; 1 -2];
y0 = [1;2];
t = [0,1];
Tol = 1e-9;

% Solving for several stepsizes 
h = [0.1; 0.01];
h_len = length(h);

 
max_err = zeros(h_len,1);


 %----------------------------------------------------------------------
 % Constant stepsize solver of third order
 %----------------------------------------------------------------------
for j = 1: h_len 
    % Constant stepsize solver
    [Ya, le, tn] = ODE23_solver(y0, t, h(j), jac, fun);

    % analytic solution
    tint = t(1):h(j):t(2);
    n = length(tint);
    Yanal = zeros(2,n);
    Yanal(:,1) = y0;
    Yanal(:, 2:end) = analLinTest(tint(2:end));

    err = zeros(1,n);
    for i = 1:n
        err(i) = norm(Yanal(:,i)-Ya(:,i));
        errmeth(i) = norm(le(:,i));
    end
    max_err(j) = max(err);
    max_errmeth(j) = max(errmeth);
end


% Checking that the order conditions are right 
ord = polyfit(log(h),log(max_err), 1);
pnum = ord(1);

% Creating convergence plot
figure; 
subplot(2,1,1);loglog(h,max_err);
hold on; axis([h(h_len)/10 h(1)*10 min(max_err)/10 max(max_err)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror advancing method for linear test problem');
hold on; 
subplot(2,1,2);loglog(h,max_errmeth,'r');
hold on; axis([h(h_len)/10 h(1)*10 min(max_err)/10 max(max_err)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror error method for linear test problem')



 %----------------------------------------------------------------------
 % Constant stepsize solver of third order
 %----------------------------------------------------------------------

 h0 = 0.1;
 
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs(fun, jac, t, y0, Tol, h0);

    % Analytic solution
    Yanal = analLinTest(tint(end));
    
    err = norm(Yanal- y(:,end));



