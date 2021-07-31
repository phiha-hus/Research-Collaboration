%Example_3_v2

% syms s
% s*E-A;
% p1 = det(s*E-A);
% % Very nice trick in MATLAB here. Computational error leads to p1 of degree
% % 3, but in fact it is of degree 1.
% disp('The matrix polynomial det(sE-A) reads')
% vpaSols = vpa(p1,6)

matrix2latex(E,'E.tex')
matrix2latex(A,'A.tex')
matrix2latex(Ad1,'Ad1.tex')
matrix2latex(Ad2,'Ad2.tex')

matrix2latex(Q,'Q.tex')
matrix2latex(Z,'Z.tex')
matrix2latex(EES,'EES.tex')
matrix2latex(AAS,'AAS.tex')
matrix2latex(AASd1,'AASd1.tex')
matrix2latex(AASd2,'AASd2.tex')

movefile('*.tex','matrices_in_Example_Ha21\')