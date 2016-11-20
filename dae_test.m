function [ yd ] = dae_test(t,y)



yd = [ -y(3);
       -4*y(2);
        y(1)+sin(t);
];



end