function [tnext, ynext, le, iflag,nfun] = onestep(f, jac,tn, yn, h, Tolit) 

[A, c, g, s]  = method;
MAXit = 100;

iflag = -1; %No iterations done yet
[m, n] = size(jac(tn,yn));


k = 0;
nfun = 0;

Y = zeros(m,s);


Y(:,1) = yn;
K = yn;

J = eye(m) - h*g*jac(tn,yn);

delY = J \ (h*g*f(tn + c(1)*h, Y(:,1)) - Y(:,1) + K);
nfun = nfun + 1; %COUNT


for i = 2:s
    sum = 0;
    Y(:,i) = Y(:,i-1);
    for j = 1:(i-1)
        sum = sum + (A(i,j)*f(tn+h*c(j),Y(:,j)));
        nfun = nfun + 1; %COUNT
    end
    K = yn + h*sum;
    k = 0; %step count
    
    while (k ~= MAXit)
        Y(:,i) = Y(:,i) + delY;
        k = k + 1;
        delY = J \ (h*g*f(tn + c(i)*h, Y(:,i)) - Y(:,i) + K);
        nfun = nfun + 1; %COUNT
       
        
        if (norm(delY) <= Tolit)
            break;
        end
        if ((k+1) == MAXit)
            tnext = tn;
            ynext = yn;
            le = 0;
            iflag = -1;
            return;
            %error('MAXit reached');
        end
    end
    Y(:,i) = Y(:,i) + delY;
    delY = 0;

    %disp(k);
end


tnext = tn + h;

ynext = Y(:,s);
le = Y(:,s) - Y(:,s-1);

if Y(:,s) ~= yn
    iflag = 1;
end


end
    
    

