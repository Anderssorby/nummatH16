function [t, y] = erk(f, t0, tend, y0, tol, h0)

n = 1;
t = [t0];
y = [y0];
h = h0;

% Our erk method
% the order of the method
p = 4;
gamma = 1/2;
s = 4;
c = [0 2*gamma 1 1];

a = [0                             0                         0                                      0;
    gamma                          gamma                     0                                      0;
    (-4*gamma+6*gamma-1)/(4*gamma) (-2*gamma+1)/(4*gamma)    gamma                                  0;
    (6*gamma-1)/(12*gamma)         (-1/(2*gamma(2*gamma-1))) (-6*gamma^2+6*gamma-1/(3*(2*gamma-1))) gamma
    ];

bt = [(-4*gamma^2+6*gamma-1)/(4*gamma) (-2*gamma+1)/(4*gamma)    gamma                                  0];

b = [(6*gamma-1)/(12*gamma)            -1/(12*gamma*(2*gamma-1)) (-6*gamma^2+6*gamma-1)/(3*(2*gamma-1)) gamma];

norm = @(x) sqrt(sum(x.^2))

% These can be recycled
asum = zeros(s);
k = zeros(s);

% local error
lerror = []


while t(n) < tend
    ssum = 0;
    for i = 1:s
        %for j = 1:s-1
        % since explicit the terms a(i,j>i)=0
        for j = 1:i
            asum(i) = asum(i) + a(i,j)*k(j);
        end
        k(i) = f(t(n) + c(i)*h, y(n) + h*asum(i));
        ssum = ssum + k(i)*b(i);
    end
    % computing the next step
    y(n+1) = y(n) + h*ssum;
    t(n+1) = t(n) + h;
    h = p*(tol/norm(lerror(n+1)))^(1/p+1)*h; 
    n = n + 1;
end
