

% Test problem
yd = @(t, y) -4*y - 4* sin(t);
z = @(t, y) [-sin(t); 4*y+4*sin(t)];

v = @(t, y) [-4*y - 4* sin(t); 
            -sin(t);
             4*y+4*sin(t)
    ];



t0 = 0;
v0 = [0; -1; 4];
h0 = 0.01;
% which gives

[Ya, Ye, tn] = ODE23_solver_mod(v0, t, h0, jac, fun);

% Checking that the order conditions are right 
ord = polyfit(log(h0),log(max_err'), 1);
pnum = ord(1)
ord = polyfit(log(h0),log(max_errmeth'), 1);
pnum = ord(1)

% Creating convergence plot
figure; 
subplot(2,1,1);loglog(h,max_err);
hold on; axis([h(h_len)/10 h(1)*10 min(max_err)/10 max(max_err)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror advancing method for the van der pol problem');
hold on; 
subplot(2,1,2);loglog(h,max_errmeth,'r');
hold on; axis([min(h)/10 max(h)*10 min(max_errmeth)/10 max(max_errmeth)*10])
xlabel('Step size');ylabel('Error');
title('Loglogerror error method for the van der pol problem')



