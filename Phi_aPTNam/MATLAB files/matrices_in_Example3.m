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
matrix2latex(Ad,'Ad.tex')

matrix2latex(hE,'hE.tex')
matrix2latex(hA,'hA.tex')
matrix2latex(hAd,'hAd.tex')

matrix2latex(P,'P.tex')
matrix2latex(K,'Abar.tex')

matrix2latex(Abar,'Abar.tex')
matrix2latex(Ad_bar,'Ad_bar.tex')

matrix2latex(hEd,'hEd.tex')
matrix2latex(hAd,'hAd.tex')

matrix2latex(H,'H.tex');
matrix2latex(H_bar,'H_bar.tex');