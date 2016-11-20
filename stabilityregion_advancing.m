%% Script plotting stability regions the error and advancing methods
%------------------------------------------------------------------------------------------
% Beautiful plot
asp = 10; b = 10;
[x, y] = meshgrid(linspace(-asp, asp), linspace(-b, b));
z = x + 1i*y;

% advancing method    
R2  = abs( -(- 6*gamma^3*z^3 + 18*gamma^2*z^3 + 18*gamma^2*z^2 - 9*gamma*z^3 - 18*gamma*z^2 - 18*gamma*z + z^3 + 3*z^2 + 6*z + 6)/(6*(gamma*z - 1)^3) ); 
    
    
%% Plots
%------------------------------------------------------------------------------------------
% Make the plot for error method
figure; clf;
contourf(x, y, R2, [1 1], 'k')
axis equal, axis([-asp asp -b b]), grid on
hold on
plot([-asp, asp], [0, 0], 'k', 'LineWidth', 1);
plot([0, 0], [-asp, asp], 'k', 'LineWidth', 1);

