% Adaptive RK-solver
% ------------------------------------------------------------------------
% Solving a system of ordinary differential equations using a RK23 method.
% ------------------------------------------------------------------------
% Input arguments:
%            fun, jac: the functions f(t, y) and Jac(t, y);
%              tn, yn: time and state variables
%                  h0: initial step size
%                 Tol: tolerance for global error.
% Output arguments:
%                   t: 
%                 le: Local error estimator.
%          iflag = 1: Iterations are successfull
%               = -1: The iterations fails. t and y are not updated
function [t, y, iflag, nfun, njac] = RKs(fun, jac, tint, y0, Tol, h0)
format long

Tolit = Tol*0.1; % Tolerance per Newton iterations
njac = 0; nfun = 0;
t = tint(1);
tend = tint(2);
it = 1;
y(:,1) = y0; %Initial conditions
h = h0;


while tend ~= t(it)
    
    % Newton iterations
    [tnext, ynext, le, iflag, njactemp, nfuntemp] = onestep(fun, jac, t(it), y(:,it), h, Tolit);
    njac = njac + njactemp; nfun = nfun + nfuntemp; % Update counters
    
    % If not converged
    if iflag == -1
        h = h*0.5;
        
    % Elseif methods not sufficiently close
    elseif norm(le) >= Tol
        h = h*0.5;
        
    % Otherwise accept 
    else
        h = 0.8 * (Tol/norm(le))^(1/4)*h;
        y(:,it+1) = ynext;
        t(it+1) = tnext;
        it = it + 1;
    end
    
    % If to long step on last iteration.
    if tnext > tend
        h = tend-t(it);
    end

end
end


