%x8 constant step size solver
     clc;
     clear all;
     close all;
     Tolit = 1e-6;
     
     format long;
     
     %Linear test function
     
     h = 0.001;
     tn = 0;
     yn = [1;2];
     jac = [-2 1; 1 -2];
     f = @(t, y) lintest(t,y);
     tend = 1;
     
     %Preallokerer stuff
     tvec = zeros(1, 100);
     yvec = zeros(2, 100);
     tvec(1) = tn;
     yvec(:,1) = yn;
     i = 1;
     while tvec(i) < tend
        [tvec(i+1), yvec(:,i+1), le, iflag] = onestep(f, jac, tvec(i), yvec(:,i), h, Tolit);
        i = i+1;
     end
     figure();
     plot(tvec, yvec);
     yAnal = zeros(2, length(tvec));
     for ndx = 1:length(tvec)
         yAnal(:, ndx) = analLinTest(tvec(ndx));
     end
     
     hold on;
     figure();
     plot(tvec, yAnal);
     
     %Van der Pol equation
     
     h = 0.001;
     tn = 0;
     yn = [0;0];
     my = 5;
     tend = 5;
     
     jac = [0,1;-my*2*yn(1)*yn(2)-1, my*(1-yn(1)^2)];
     
     v = @(t,y) vanDerPol(t,y, my);
     
          %Preallokerer stuff
     tvec = zeros(1, 100);
     yvec = zeros(2, 100);
     tvec(1) = tn;
     yvec(:,1) = yn;
     i = 1;
     while tvec(i) < tend
        [tvec(i+1), yvec(:,i+1), le, iflag] = onestep(f, jac, tvec(i), yvec(:,i), h, Tolit);
        i = i+1;
     end
     figure();
     plot(tvec,yvec(1,:));
     