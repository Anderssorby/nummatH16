%% Script plotting stability regions the error and advancing methods
%------------------------------------------------------------------------------------------
% Beautiful plot
asp = 10; b = 10;
[x, y] = meshgrid(linspace(-asp, asp), linspace(-b, b));
z = x + 1i*y;
g = 0.4395;

% error method
R1  = abs(((z2.^2*(2*g^2-4*g+1)+z2.*(-4*g+2)+2)./(2*(1-g*z2).^2)));



%% Plots
%------------------------------------------------------------------------------------------
% Make the plot for error method
figure; clf;
contourf(x, y, R1, [1 1], 'k')
axis equal, axis([-asp asp -b b]), grid on
hold on
plot([-asp, asp], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-asp, asp], 'k', 'LineWidth', 1);
title('Stabilityregion for error method');
xlabel('real axis');
ylabel('imaginary axis');

