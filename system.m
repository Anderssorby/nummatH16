syms gamma
s = 4;
c = [0 2*gamma 1 1];

a = [0                             0                         0                                      0;
    gamma                          gamma                     0                                      0;
    (-4*gamma^2 + 6*gamma - 1)/(4*gamma) (-2*gamma + 1)/(4*gamma)    gamma                                  0;
    (6 * gamma - 1)/(12 * gamma)         (-1/(12*gamma*(2*gamma - 1))) (-6*gamma^2 + 6*gamma - 1)/(3*(2*gamma - 1)) gamma
    ];

bt = [(-4*gamma^2+6*gamma-1)/(4*gamma) (-2*gamma+1)/(4*gamma)    gamma                                  0];

b = [(6*gamma-1)/(12*gamma)            -1/(12*gamma*(2*gamma - 1)) (-6*gamma^2 + 6*gamma - 1)/(3*(2*gamma - 1)) gamma];


