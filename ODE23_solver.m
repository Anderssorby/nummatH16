% Constant stepsize RK-solver
% Solving a system of ordinary differential equations using an embedded
% pair of implicit methods. One second order and one third order method.

function [Ya, Ye, tn] = ODE23_solver(y0, t, h , jac, fun)
format long

Tolit = 1e-9; % Tolerance for Newton iterations
tn = t(1);
tend = t(2);
iterations = ceil((tend-tn) / h); 
Ye = zeros(length(y0),iterations+1); % preallocation error
Ya = zeros(length(y0),iterations+1); % preallocation advancing meth
%Initial conditions
Ya(:,1) = y0;

for i = 1:iterations
    % Newton iterations
    [tnext, ynext, le, iflag] = onestep(fun, jac, tn, Ya(:,i), h, Tolit);
    if iflag == -1
        error('Too many newton iterations');
    end
    Ye(:,i+1) = le;
    Ya(:,i+1) = ynext;
    tn = tnext;
end
end


