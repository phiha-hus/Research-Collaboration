function f=circuit_1_d1(t,i,di)

L = 1/3;
C = 1.5;
R = 2;

f = [ -1 0 0 0
      0 0 0 -1/C
      0 0 1 -R
      1 1 1 0];
