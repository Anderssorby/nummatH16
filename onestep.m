%% Newton fixed point iteration 
% ------------------------------------------------------------------------
% Do one step with an implicit RK23 method. The formulas referred to is
% given in the project description.
% ------------------------------------------------------------------------
% Input arguments:
%              f, jac: the functions f(t, y) and Jac(t, y);
%              tn, yn: time and state variables
%                   h: step size
%               Tolit: tolerance for the Newton iterations.
% Output arguments:
%       tnext, ynext: time and state variables after one step
%                 le: Local error estimator.
%          iflag = 1: Iterations are successfull
%               = -1: The iterations fails. t and y are not updated
%
% ------------------------------------------------------------------------
%% [tnext, ynext, le, iflag, njac, nfun] = onestep(f, jac, tn, yn, h, Tolit)
% ------------------------------------------------------------------------
function [tnext, ynext, le, iflag, njac, nfun] = onestep(f, jac, tn, yn, h, Tolit)
format long
    
   [A, c, g, s] = method; % Embedded third order RK 
   maxstep = 100;
   len_yn = length(yn);
   Y = zeros(len_yn,s); % solution
   Y(:, 1) = yn;  % initial value
   Idyy = eye(size(jac(tn,yn)));
   njac = 1; nfun = 0; % Counters
   iflag = 1;
   
   % For every k
   for i = 2:s
        ksum = 0; %preallocate summation
        Y(:,i) = Y(:,i-1); % Update initial New
        
        % Summation (3b)
        for j = 1:(i-1)
            ksum = ksum + A(i,j)*f(tn + c(j)*h, Y(:,j));
            nfun = nfun +1; % Count
        end
        
        k = yn + h*ksum;
        it = 0; % Step count
        
        % Newton fixed point iteration
        while it <= maxstep
            
            lhs = Idyy - h*g*jac(tn+c(i)*h,Y(:,i)); % jacobian update
            
            % Newton iteration (4) 
            del_Y = lhs\(h*g*f(tn+c(i)*h,Y(:,i))-Y(:,i)+k); 
            Y(:,i) = Y(:,i) + del_Y; 
            it = it + 1; njac = njac +1;nfun = nfun +1; % Counters
            
            if norm(del_Y) <= Tolit % If converged
                break;
            end
        end
        
        % if too many iterations
        if it > maxstep
            % Return variables.
            ynext = yn;
            le = 0;
            tnext = tn;
            iflag = -1;
            return;
        end
   end
    
    % Return variables.
    ynext = Y(:, s); % update (3c)
    le = Y(:, s) - Y(:,s-1); % Update (3d)
    tnext = tn + h; 
end