a = 10; 
b = 10;
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
z = x + 1i*y;

% Stability function.
R_ad_abs= abs((z.^3*(-12*gamma^4+42*gamma^3-36*gamma^2+11*gamma-1)+z.^2*...
    (36*gamma^3-54*gamma^2+24*gamma-3)+z.*(-36*gamma^2+30*gamma-6)...
    +12*gamma-6)./(6*(2*gamma-1)*(1-gamma*z).^3));

% Make the plot.
figure(1)
contourf(x, y, R_ad_abs, [1 1], 'k')
axis equal, axis([-a a -b b]);
grid on; hold on;
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);