clear all; close all

a = 10; 
b = 10;
g = 0.4358665215; %gamma found in stabilityFunctions.m
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
z = x + 1j*y;

% Stability function for the advancing method
R = abs(((z.^3*(-12*g^4+42*g^3-36*g^2+11*g-1)+z.^2*(36*g^3-54*g^2+24*g-3)+...
    z.*(-36*g^2+30*g-6)+12*g-6)) ./(6*(2*g-1)*(1-g*z).^3));

c = 100;
d = 100;
[x2, y2] = meshgrid(linspace(-c, c), linspace(-d, d));
z2 = x2 + 1j*y2;
g = 0.4358665215; %gamma found in stabilityFunctions.m

%stability function for the error estimating method
R_hat = abs(((z2.^2*(2*g^2-4*g+1)+z2.*(-4*g+2)+2)./(2*(1-g*z2).^2)));


% Plotting the stability regions
figure(1)
subplot(1,2,1)
contourf(x, y, R, [1 1], 'k')
axis equal, axis([-a a -b b]);
grid on; hold on;
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);
xlabel('Im(z)')
ylabel('Re(z)')
title('Stability region for the advancing method')
subplot(1,2,2)
contourf(x2, y2, R_hat, [1 1], 'k')
axis equal, axis([-c c -d d]);
grid on; hold on;
plot([-c, c], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-c, c], 'k', 'LineWidth', 1);
xlabel('Im(z)')
ylabel('Re(z)')
title('Stability region for the error estimating method')
