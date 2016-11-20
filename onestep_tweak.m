%% Newton fixed point iteration 
% ------------------------------------------------------------------------
% Do one step with an implicit RK23 method on a DEA. The formulas referred to is
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
%% [tnext, ynext, le, iflag, njac, nfun] = onestep(f, jac, M, tn, yn, h, Tolit)
% ------------------------------------------------------------------------
function [tnext, ynext, le, iflag, njac, nfun] = onestep_tweak(f, jac, M, tn, vn, h, Tolit)
format long
    
   [A, c, g, s] = method; % Embedded third order RK 
   maxstep = 100;
   len_yn = length(vn);
   V = zeros(len_yn,s); % solution
   V(:, 1) = vn;  % initial value
   njac = 0; nfun = 0; % Counters
   iflag = 1;
   
   % For every k
   for i = 2:s
        ksum = 0; %preallocate summation
        V(:,i) = V(:,i-1); % Update initial New
        
        % Summation (3b)
        for j = 1:(i-1)
            ksum = ksum + A(i,j)*f(tn + c(j)*h, V(:,j));
            nfun = nfun +1; % Count
        end
        
        k = M*vn + h*ksum;
        it = 0; % Step count
        
        % Newton fixed point iteration
        while it <= maxstep
            
            lhs = M - h*g*jac(tn+c(i)*h,V(:,i)); % jacobian update
            
            % Newton iteration (4) 
            del_Y = lhs\(h*g*f(tn+c(i)*h,V(:,i))-M*V(:,i)+k); 
            V(:,i) = V(:,i) + del_Y; 
            it = it + 1; njac = njac +1;nfun = nfun +1; % Counters
            
            if norm(del_Y) <= Tolit % If converged
                break;
            end
        end
        
        % if too many iterations
        if it > maxstep
            % Return variables.
            ynext = vn;
            le = 0;
            tnext = tn;
            iflag = -1;
            return;
        end
   end
    
    % Return variables.
    ynext = V(:, s); % update (3c)
    le = V(:, s) - V(:,s-1); % Update (3d)
    tnext = tn + h; 
end