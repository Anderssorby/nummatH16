% Constant stepsize RK-solver
% Solving a system of ordinary differential equations using an embedded
% pair of implicit methods. One second order and one third order method.
% - Modified for linear systems: Mv' = F(t,v)

function [Ya, Ye, tn] = ODE23_solver_mod(v0, t, h , jac_g, g, f)
format long

Tolit = 1e-7; % Tolerance for Newton iterations

% solve for z0 and y0
% Newton iterations

v = v0;

it = 0; % Step count
lhs = Idyy - h*g*jac(tn+c(i)*h,Y(:,i)); % jacobian calc
temp = Y(:,i)+Tolit*10; % 

% Newton fixed point iteration
while max(abs(Y(:,i)- temp)) > Tolit % Jac change
            % if too many steps
     if it >= maxit
        iflag = -1;
        tnext = tn;
        ynext = yn;
        le = yn;
     end

    temp = Y(:,i); % previous iteration
    Y(:,i) = temp + lhs\(h*g*fun(tn+c(i)*h,temp)-temp+k); % func calc
    it = it + 1; nfun = nfun + 1; njac = njac + 1; % Counters
   
    lhs = Idyy - h*g*jac(tn+c(i)*h,temp); % jacobian calc
end



tn = t(1);
tend = t(2);
iterations = ceil((tend-tn) / h); 
Ye = zeros(length(y0),iterations+1); % preallocation error
Ya = zeros(length(y0),iterations+1); % preallocation advancing meth
%Initial conditions
Ya(:,1) = y0;

for i = 1:iterations
    % Newton iterations
    [tnext, ynext, le, iflag] = onestep_adapt_mod(fun, jac, tn, Ya(:,i), h, Tolit);
    if iflag == -1
        error('Too many newton iterations');
    end
    Ye(:,i+1) = le;
    Ya(:,i+1) = ynext;
    tn = tnext;
end
end


