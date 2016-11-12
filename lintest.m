function yut = lintest(t,y)
    yut = zeros(2,1);
    yut(1) = -2*y(1) + y(2)+ t;
    yut(2) = y(1) -2*y(2) + t + 3;
end