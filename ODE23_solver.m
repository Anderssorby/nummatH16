% Constant stepsize RK-solver
% Solving a system of ordinary differential equations using an embedded
% pair of implicit methods. One second order and one third order method.

function [Y4, Ye, Y3, tn] = ODE23_solver(y0, t, h , jac, fun)
format long

Tolit = 1e-9; % Tolerance for Newton iterations
tn = t(1);
tend = t(2);
iterations = ceil((tend-tn) / h); 
Ye = zeros(length(y0),iterations+1); % preallocation error
Y4 = zeros(length(y0),iterations+1); % preallocation advancing meth
Y3 = zeros(length(y0),iterations+1); % preallocation lower order meth
%Initial conditions
Y4(:,1) = y0;

for i = 1:iterations
    % Newton iterations
    [tn, Y4(:,i+1), Y3(:,i+1), Ye(:,i+1), iflag] = onestep(fun, jac, tn, Y4(:,i), h, Tolit);
    if iflag == -1
        error('Too many newton iterations');
    end
end
end


