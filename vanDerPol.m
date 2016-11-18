function yder = vanDerPol(t,y, my)
    yder = zeros(2,1);
    yder(1) = y(2);
    yder(2) = my*(1-y(1)^2)*y(2)-y(1);
end