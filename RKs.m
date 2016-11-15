% Adaptive RK-solver
% Solving a system of ordinary differential equations using an embedded
% pair of implicit methods. One second order and one third order method.

function [t, y, iflag, nfun, njac] = RKs(fun, jac, tint, y0, Tol, h0)
format long

Tolit = Tol*0.1; % Tolerance for Newton iterations
njac = 0; nfun = 0;
t = tint(1);
tend = tint(2);
it = 1;
%Initial conditions
y(:,1) = y0;
iflag = -1;


while abs(t(it) - tend) > Tolit
    % Newton iterations
    [tnext, ynext, le, iflag, njactemp, nfuntemp] = onestep_adapt(fun, jac, t(it), y(:,it), h0, Tolit);
    njac = njac + njactemp; nfun = nfun + nfuntemp; % Update counters
    while iflag == -1
        h = h0*0.5; 
        [tnext, ynext, le, iflag, njactemp, nfuntemp] = onestep_adapt(fun, jac, t(it), y(:,it), h, Tolit);
        njac = njac + njactemp; nfun = nfun + nfuntemp; % Update counters
    end
    if tnext > tend
        [tnext, ynext, le, iflag, njactemp, nfuntemp] = onestep_adapt(fun, jac, t(it), y(:,it), tend-t(it), Tolit);
        njac = njac + njactemp; nfun = nfun + nfuntemp; % Update counters
    end
    y(:,it+1) = ynext;
    t(it+1) = tnext;
    it = it + 1;
end
njac
nfun
end


