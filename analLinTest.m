function y = analLinTest(t)
    %Gives analytic solution to linear test problem
    len = length(t);
    y = zeros(2,len);
    
    y(1,:) = exp(-t) + t;
    y(2,:) = exp(-t) + t + 1;
end