% Adaptive RK-solver
% Solving a system of ordinary differential equations using an embedded
% pair of implicit methods. One second order and one third order method.

function [t, y, iflag, nfun, njac] = RKs(fun, jac, tint, y0, Tol, h0)
format long

Tolit = Tol*0.05; % Tolerance for Newton iterations
njac = 0; nfun = 0;
t = tint(1);
tend = tint(2);
it = 1;
%Initial conditions
y(:,1) = y0;
iflag = -1;

%Safety constant for h-reduction
safeFactor = 0.69;
%An attempt at finishing at the right place
finishDist = Tolit;


while (t(it) - tend) < 0
    % Newton iterations
    
    h = h0;
    
    [tnext, ynext, le, iflag, njactemp, nfuntemp] = onestep_adapt(fun, jac, t(it), y(:,it), h, Tolit);
    njac = njac + njactemp; nfun = nfun + nfuntemp; % Update counters
    while iflag == -1 || (tnext - tend) > finishDist
        
        %masse lok for å finne smart måte å ikke oversteppe
        if (tnext - tend) > finishDist 
            h = tend-t(it);
        else
            h = safeFactor*Tol/norm(le)^(1/4); 
        end
        
        
        [tnext, ynext, le, iflag, njactemp, nfuntemp] = onestep_adapt(fun, jac, t(it), y(:,it), h, Tolit);
        njac = njac + njactemp; nfun = nfun + nfuntemp; % Update counters
        
    end
    t(it+1) = tnext;
    y(:,it+1) = ynext;
    it = it + 1;
end
njac
nfun
end


