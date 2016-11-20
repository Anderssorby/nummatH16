%% Main linear test problem
% ------------------------------------------------------------------------
% Solving linear test problem with constant step size 

clear all
close all

%% variable definition
fun = @(t, y) lintest(t,y);
jac = @(t, y) lintest_jac(t,y);
y0 = [1;2]; % Initial values
tint = [0,1]; % Time interval
Tol = 1e-7; % Global tolerance
h = [0.1; 0.01; 0.001]; % Step sizes
h_len = length(h);

 
advmeth = zeros(h_len,1);
errmeth = zeros(h_len,1);



%% Solver 
for j = 1: h_len 
    % Constant stepsize solver
    [y, le, t, njac, nfun] = ODE23_solver(y0, tint, h(j) , jac, fun);
    
    % Analytic solution
    yanal = lintest_anal(tint(2));


    advmeth(j) = max(y(:,end) - yanal(:)); % Advancing method
    errmeth(j) = max(max(le)); % Error method 
end


%% Postprocessing
% Error conditions
ord = polyfit(log(h),log(advmeth), 1);
pnum_adv = ord(1)
ord = polyfit(log(h),log(errmeth), 1);
pnum_err = ord(1)

% Second order line
line_2 = zeros(2,1);
line_2(1) = h(1)^2/10;
line_2(2) = h(end)^2/10;
hvec_2 = [h(1);h(end)];

% Third order line
line_3 = zeros(2,1);
line_3(1) = h(1)^3/100;
line_3(2) = h(end)^3/100;
hvec_3 = [h(1);h(end)];


% ------------------------------------------------------------------------
% Convergence plots
% ------------------------------------------------------------------------
% Error method
figure(); grid on;
subplot(2,1,1);loglog(h,errmeth,'b'); 
hold on; loglog(hvec_2, line_2,'r');
%hold on; axis([h(h_len)/10 h(1)*10 min(errmeth)/10 max(errmeth)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror error method for linear test problem')

% Advancing method
subplot(2,1,2);loglog(h,advmeth, 'b'); 
hold on; loglog(hvec_3, line_3, 'r');
%hold on; axis([h(h_len)/10 h(1)*10 min(advmeth)/10 max(advmeth)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror advancing method for linear test problem')

