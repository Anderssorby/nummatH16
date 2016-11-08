% Plot the stabilty domain for a given stability function.
% Domain: [−a, a, −b, b]
a = 4; b = 4;
[x, y] = meshgrid(linspace(-a, a), linspace(-b, b));
z = x + i*y;
% Stability function.
R = abs(1 + z + z.^2/2);
% Make the plot.
contourf(x, y, R, [1 1], 'k')
axis equal, axis([-a a -b b]), grid on
hold on
plot([-a, a], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-a, a], 'k', 'LineWidth', 1);