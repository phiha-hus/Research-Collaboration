% Anh nho Muoi giup anh code thu xem sao nhe

E = blkdiag(eye(3),zeros(3,3))
     
     
Abar = [-0.0050    0.0005   -0.0005
     -0.0083    0.0008   -0.0008
     0          0         0];

A = blkdiag(Abar,eye(3))

Ad_bar = [0.0000    0.0007    0.0067
        0.0001    0.0011    0.0112
         0         0         0];
     
K =[0.1523    0.1530    1.5551
    1.5781    0.2716    2.5919
         0         0         0];
     
Ad = [Ad_bar Ad_bar;K K]     