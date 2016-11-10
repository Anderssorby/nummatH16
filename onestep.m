function [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit)
% [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit)
% Do one step with an implicit RK−method method.
% Input arguments:
% f, jac: the functions f(t, y) and Jac(t, y);
% tn, yn: time and state variables
% h: step size
% Tolit: tolerance for the Newton iterations.
% Output arguments:
% tnext, ynext: time and state variables after one step
% le: Local error estimator.
% iflag = 1: Iterations are successfull
% = −1: The iterations fails. t and y are not updated
[A, g, c, s] = method()
Y = zeros(s)
Y(1) = yn;
for i = 2:s
    ksum = 0;
    for j = 1:s
        ksum = ksum + A(i,j)*f(tn + h*c(i), Y(j));
    end
    
    cons = yn + h*ksum; 
    
    % Newton iterations
    N = cons; % better start value?
    it = 0;
    while it < Tolit % count eveluations of jac
        N = N - jac(tn + h*c(i), N)\(cons + h*f(tn + h*c(i), N)); 
        it = it + 1;
    end
    Y(i) = N;
end
%accept?
iflag = 1;
le = Y(4) - Y(3);
ynext = Y(4);

tnext = tn + h

end
