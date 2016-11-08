
% Our erk method
% the order of the method
p = 4; %?
syms gamma
s = 4;
c = [0 2*gamma 1 1];

a = [0                             0                         0                                      0;
    gamma                          gamma                     0                                      0;
    (-4*gamma^2 + 6*gamma - 1)/(4*gamma) (-2*gamma + 1)/(4*gamma)    gamma                                  0;
    (6 * gamma - 1)/(12 * gamma)         (-1/(12*gamma*(2*gamma - 1))) (-6*gamma^2 + 6*gamma - 1)/(3*(2*gamma - 1)) gamma
    ];

bt = [(-4*gamma^2+6*gamma-1)/(4*gamma) (-2*gamma+1)/(4*gamma)    gamma                                  0];

b = [(6*gamma-1)/(12*gamma)            -1/(12*gamma*(2*gamma - 1)) (-6*gamma^2 + 6*gamma - 1)/(3*(2*gamma - 1)) gamma];

norm = @(x) sqrt(sum(x.^2));


% p = 1 a single node



% p = 3
%fprintf('p = 3, b = %s\n', b)
phi = 0;
for i=1:s
    phi = phi + b(i)*c(i)*c(i);
end
simplify(phi)

phi = 0;
for i=1:s
    for j=1:s
        phi = phi + b(i)*a(i, j)*c(j);
    end
end
simplify(phi)




%fprintf('p = 3, bt = %s\n', bt)
phi = 0;
for i=1:s
    phi = phi + bt(i)*c(i)*c(i);
end
fprintf('computing for bt')
simplify(phi)

phi = 0;
for i=1:s
    for j=1:s
        phi = phi + bt(i)*a(j, i)*c(j);
    end
end

simplify(phi)

