% RK-solver
% Solving a system of ordinary differential equations using an embedded
% pair of implicit methods. One second order and one third order method.



h = 0.1; % constant step
Tolit = 1e-6;
jac = [-2 1; 1 -2];
tn = 0;
yn = [1;2];
tend = 10;
iterations = (tend-tn) / h; 
Ye = zeros(length(yn),iterations+1); % preallocation error
Ya = zeros(length(yn),iterations+1); % preallocation advancing meth
Yanal = zeros(length(yn),iterations+1); % preallocation analytic solution
Ya(:,1) = yn;
Yanal(:,1) = yn;

% Problem
f = @(t, y) lintest(t,y);

for i = 1:iterations
    [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit);
    Ye(:,i+1) = le;
    Ya(:,i+1) = ynext;
    tn = tn + tnext;
    Yanal(:,i) = analLinTest(tn); % analyticv solution
end



