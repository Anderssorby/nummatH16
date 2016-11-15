function [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit)
format long
    % [tnext, ynext, le, iflag] = onestep(f, jac, tn, yn, h, Tolit)
    % Do one step with an implicit RK?method.
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
    
    % RK third order method
   [A, c, g, s] = method;
    
   maxstep = 100;
   len_yn = length(yn);
   Y = zeros(len_yn,4); % solution
   Y(:, 1) = yn;  % initial value
   Idyy = eye(size(jac(tn,yn)));
   iflag = 1;
  
   % For every k
   for i = 2:s
        summa = zeros(len_yn,1); %preallocate summation
        
        % For sum in RK
        for j = 1:i-1
            summa = summa + A(i,j)*f(tn + c(j)*h, Y(:,j));
        end
        
        k = yn + h*summa;
        it = 0; % Step count
        
        temp = Y(:,i)+Tolit*10; % 
        
        % Newton fixed point iteration
        while max(abs(Y(:,i)- temp)) > Tolit && it < maxstep % Jac change
            
            temp = Y(:,i); % previous iteration
            lhs = Idyy - h*g*jac(tn+c(i)*h,temp); % jacobian calc
            Y(:,i) = temp + lhs\(h*g*f(tn+c(i)*h,temp)-temp+k);
            it = it + 1;
        end
        % if too many steps
        if it > maxstep
            iflag = -1;
            return;
        end
   end
    
    % Return variables.
    ynext = Y(:, 4);
    le = Y(:, 4) - Y(:,3);
    tnext = tn + h;
end