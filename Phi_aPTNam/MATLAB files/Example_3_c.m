%Example_3_v2
% The command Example_3_v2 is ALREADY IN matrices_in_Example3 
matrices_in_Example3

% Example_3_b >>> Linear Programming does not give 
% us trustable results.

[Up,Sp,Vp] = svd(P);
Sp

A1 = Up' * Abar * Vp

% H1 = Up' * H * Up has the same stability property
... as H
    
H1 = blkdiag(A1(1,1)/Sp(1,1),diag(-rand(2,1))) ;

H = Up * H1 * Up'
eig(H)  %-0.5499    -0.3264    -0.1304

n = 3;
I = eye(n);

H_bar = [Ad_bar+H Ad_bar; K K-I] ;
eig(H_bar)
% -1.2832    -0.1221    -0.3264    -0.3001    -0.5499    -1.0000

matrix2latex(H,'H.tex');
matrix2latex(H_bar,'H_bar.tex');