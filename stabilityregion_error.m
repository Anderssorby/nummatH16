%% Script plotting stability regions the error and advancing methods
%------------------------------------------------------------------------------------------
% Beautiful plot
asp = 10; b = 10;
[x, y] = meshgrid(linspace(-asp, asp), linspace(-b, b));
z = x + 1i*y;


% error method
R1  = abs( (2*gamma^2*z^2 - 4*gamma*z^2 - 4*gamma*z + z^2 + 2*z + 2)/(2*(gamma*z - 1)^2) ); 



%% Plots
%------------------------------------------------------------------------------------------
% Make the plot for error method
figure; clf;
contourf(x, y, R1, [1 1], 'k')
axis equal, axis([-asp asp -b b]), grid on
hold on
plot([-asp, asp], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-asp, asp], 'k', 'LineWidth', 1);

