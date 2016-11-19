function y = lintest_anal(t)
    %Gives analytic solution to linear test problem
    y = zeros(2,length(t));
    y(1, :) = exp(-t) + t;
    y(2, :) = exp(-t) + t +1;
end