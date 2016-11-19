function [tnext, ynext, le, iflag, njac, nfun] = onestep_adapt_mod(fun, jac, tn, yn, h, Tolit)
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
    
   njac = 1; % Counter evaluations of jac
   nfun = 0; % Counter evaluations of f
   maxit = 10; % Maximum number of fixed point iterations
   len_yn = length(yn);
   Y = zeros(len_yn,4); % solution
   Y(:, 1) = yn;  % initial value
   Idyy = eye(size(jac(tn,yn))); % jac calc
   iflag = 1;
  
   % For every k
   for i = 2:s
        %summa = zeros(len_yn,1); %preallocate summation
        summa = 0;        

        % For sum in RK
        for j = 1:i-1
            summa = summa + A(i,j)*fun(tn + c(j)*h, Y(:,j)); % fun calc
            nfun = nfun + 1;
        end
        
        k = yn + h*summa;
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
                return;
             end

            temp = Y(:,i); % previous iteration
            Y(:,i) = temp + lhs\(h*g*fun(tn+c(i)*h,temp)-temp+k); % func calc
            it = it + 1; nfun = nfun + 1; njac = njac + 1; % Counters
           
            lhs = Idyy - h*g*jac(tn+c(i)*h,temp); % jacobian calc
        end

   end
    % Return variables.
    ynext = Y(:, 4);
    le = Y(:, 4) - Y(:,3);
    tnext = tn + h;
    
    if max(abs(le)) > Tolit
        iflag = -1;
    end
end