format long
rksystem()

%syms la h z

%z = la*h

%Rhla = (1+a(3,1)*h*la)/(1-gamma*h*la) + (a(3,2)*h*la + a(3,3)*h*la*(1+a(2,1)*h*la))/((1-gamma*h*la)^2) + (a(3,3)*a(2,2)*(h*la)^2)/((1-gamma*h*la)^3)

%R(z) = (1 + a(3,1) * z)/(1 - gamma * z) + ...
%    (a(3,2) * z + a(3,3) * z * (1 + a(2,1) * z))/((1-gamma*z)^2) + ...
%    (a(3,3)*a(2,2)*z^2)/((1-gamma*z)^3);


gamma = 0.29289;
%R = subs(R)

asp = 2; b = 2;
n = 50;
[x, y] = meshgrid(linspace(-asp, asp, n), linspace(-b, b, n));
z = x + i*y;
% Stability function.
a = double(subs(a));
R = abs( (1 + gamma*z)/(1 - gamma*z) + ...
    (a(3,2)*z*(1 + a(2,1)*z) + a(4,3)*z*(1 + a(3,1)*z))/((1 - gamma*z)^2) + ...
    a(4,3)*a(3,2)*z^2*(1 + a(2,1)*z)/((1 - gamma*z)^2) );
Rtest = abs(1 + z);
% Make the plot.
figure; clf;
contourf(x, y, R, [1 1], 'k')
axis equal, axis([-asp asp -b b]), grid on
hold on
plot([-asp, asp], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-asp, asp], 'k', 'LineWidth', 1);
