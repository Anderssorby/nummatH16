% Main Circuit differentiator
% ------------------------------------------------------------------------
% Solving the circuit differentiator DAE with adaptive step size

clear all
close all

%% Variable definitions

% Physical variables 
C = 1e-12; % [F]
V = @(t) voltage(t); % [V]
R = 1e4; % [ohm]
A = 300;

% System properties
fun = @(t, y) dae_circuit(t, y, V, R, A,C);
jac = @(t, y) dae_circuit_jac(t, y, R, A,C);
fun_2 = @(t, y) dae_circuit_2(t, y, V, R, A);
jac_2 = @(t, y) dae_circuit_jac_2(t, y, R, A);
M = [C -C 0 0 0; -C C 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0];
y0 = [0]; % Initial conditions. 
y0_2 = [0;0;0;0;0];
tint = [0,1e-8]; % Time interval
Tol = 1e-10;
len_Tol = length(Tol); 
h0 = 1e-9; % Initial first step 


%% Solver
for i = 1:len_Tol
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs(fun, jac, tint, y0, Tol(i), h0);
    
    % Work
    werk(i) = nfun + 3*njac;
end
      
z2 = V(t)-y(1,:);
z3 = -A*(V(t)-y(1,:));
z4 = -(1+A)*z2/R;
z5 = -(1+A)*(V(t)-y)/R;

figure();
hold on;
subplot(2,3,1);plot(t,y(1,:))
title('V(t)');
hold on;
subplot(2,3,2);plot(t,z2(:))
title('u_2');
hold on;
subplot(2,3,3);plot(t,z3(:))
title('u_3');
hold on;
subplot(2,3,4);plot(t,z4(:))
title('I_a');
hold on;
subplot(2,3,5);plot(t,z5(:)) 
title('I_Q');
hold on;
subplot(2,3,6);
text(0,0.65,'ODE form of problem'); axis off;
%% Solver
for i = 1:len_Tol
    % Adaptive stepsize solver
    [t,y,iflag,nfun,njac] = RKs_tweak(fun_2, jac_2, M, tint, y0_2, Tol(i), h0);
    
    % Work
    werk(i) = nfun + 3*njac;
end

figure();
hold on;
subplot(2,3,1);plot(t,V(t), 'm');
hold on;
title('V(t)');
a = subplot(2,3,2);plot(t,y(2,:),'m')
hold on;
title('u_2');
subplot(2,3,3);plot(t,y(3,:),'m')
hold on;
title('u_3');
subplot(2,3,4);plot(t,y(4,:),'m')
hold on;
title('I_a');
subplot(2,3,5);plot(t,y(5,:),'m')
hold on;
title('I_Q');
subplot(2,3,6);
text(0,0.65, 'DAE form of problem'); axis off;