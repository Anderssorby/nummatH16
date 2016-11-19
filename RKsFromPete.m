function [t, y, iflag, nfun, njac, hvec] = RKs(f, jac, t0, tend, y0, Tol, h0)

Tolit = 0.1 * Tol;

[m,n] = size(jac(t0,y0));

nfun = 0;
njac = 0;


%Y = zeros(m,n);
%Yhat = zeros(m,n);

P = 0.69;
p = 3;

ind = 0;

Y(:,1) = y0;
t(1) = t0;

i = 1;


ynext = y0;
tnext = t0;
h = h0;

%[tnext, ynext, le, iflag] = onestepfuck(f, jac,t0, y0, h0, Tolit);


while abs(tnext-tend) > 1e-14
    
    if h+tnext > tend
        h = tend-tnext;
    end
    

    [ttemp, ytemp, le, iflag,ntemp] = onestepFromPete(f, jac,tnext, ynext, h, Tolit);
    nfun = nfun + ntemp; %teller hver gang den kjører
    njac = njac + 1;
    
    if iflag == -1
       h = h/2;
       disp('Ingen konvergens')
       disp('i er: ')
       disp(i);
       disp('------')
    elseif norm(le) <= Tol
       i = i + 1;
       ynext = ytemp;
       tnext = ttemp;
       
       Y(:,i) = ynext;
       t(i) = tnext;
       hvec(i) = h;
       
       h = h0;
   
    else
       h = P * (Tol/norm(le))^(1/(p+1))*h;
       
       disp('redusering notert');
       disp(h);
       disp('i er: ');
       disp(i);
       disp('------');
       
    end
    

end

% 
% 
% 
% disp('antall steg: ');
% disp(i);
% 
% disp(length(t));
% disp(length(Y(1,:)));

%{
hold on
for i = 1:m
    plot(t,Y(i,:),'Color',colormaker(i));
    legendInfo{i} = ['Y = ' num2str(i)];
end
legend(legendInfo);
grid on;
%}

y = Y;

end