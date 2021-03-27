%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% definition of matrix for truck example
%
% CALL  : M=example_truck(t,T,h,y,yi,K,a,sel)
%
%
% INPUT : t   - current time t\in[t0,t0+a]
%         T   - time steps T=[t0,t1,..,t]
%         h   - current stepsize
%         y   - vector of current values
%         yi  - vector of previous values
%         K   - current order
%         a   - vector of bdf coefficients
%         sel - selection
%                 sel=1 - dimension of [q,r,w,s,lambda]
%                 sel=2 - evaluation of F
%                 sel=3 - evaluation of DF
%
% OUTPUT: M   - matrix or vector 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function f = truck(t,y,dy,ddy)
   
p=y(1:11);
lambda=10e-6*y(12);
pd=dy(1:11);
%lambdad=y(24);

cosp4= cos(p(4));
sinp4= sin(p(4));
cosp6= cos(p(6));
sinp6= sin(p(6));
cosp8= cos(p(8));
sinp8= sin(p(8));
cosp11= cos(p(11));
sinp11=sin(p(11));

g=9.81;

m1=1450;
m2=600;
m3=1980;
m4=1355;
m5=1000;
m6=100;
m7=11515;
l3=10367;
l4=432;
l5=948;
l7=33000;

k10=4400000;
k20=2200000;
k34=2643833;
k43=779735;
k35=135707;
k53=135707;
k56=1000;
k37=900000;

d10=600;
d20=300;
d34=618;
d43=182;
d35=12218;
d53=12218;
d56=447;
d37=38500;

a13=2.0625;
a23=2.4375;
a35=1.9375;
a34=1.8125;
a43=3.0375;
a53=3.6375;
a37=0.9825;
az1=3.0625;
az2=0.1500;

b34=0.279;
b43=0.946;
c35=0.800;
c53=0.900;
c56=0.200;
e37=2.435;
ez1=1.610;
ez2=0.750;
heq=0.900;


Feq10 =        0.1226542774000000E+06;           
Feq20 =        0.5392572260000000E+05;            
Feq34 =        0.1026510391836735E+05;        
Feq43 =        0.3027446081632654E+04;           
Feq35 =        0.5597470588235294E+04;        
Feq53 =        0.5193529411764707E+04;          
Feq56 =        0.9810000000000000E+03;           
Feq37 =       -0.6300385509270705E+06;   
Feq13 =       -0.1084297774000000E+06;     
Feq23 =       -0.4803972260000000E+05;  


P=1.0e+5;
A=0.0562;
s=39.1328;
k=1.4;
d13=21593;
d23=38537;
lnom13=0.160;
lnom23=0.160;
l=0.360;

t1=0.15;
t2=0.0;

x10=p(1)-excitation(t-t1);
x20=p(2)-excitation(t-t2);
x13=p(3)+a13*sin(p(4))-p(1);
x23=p(3)-a23*sin(p(4))-p(2);
x34=p(5)+b34*sin(p(6))-p(3)+a34*sin(p(4));
x43=p(5)-b43*sin(p(6))-p(3)+a43*sin(p(4));
x35=p(7)+c35*sin(p(8))-p(3)+a35*sin(p(4));
x53=p(7)-c53*sin(p(8))-p(3)+a53*sin(p(4));
x56=p(9)-p(7)+c56*sin(p(8));
x37=p(10)-e37*sin(p(11))-ez2*cos(p(11))-p(3)+a37*sin(p(4));
xp10=pd(1)-excitation_p(t-t1);
xp20=pd(2)-excitation_p(t-t2);
xp13=pd(3)+a13*pd(4)*cos(p(4))-pd(1);
xp23=pd(3)-a23*pd(4)*cos(p(4))-pd(2);
xp34=pd(5)+b34*pd(6)*cos(p(6))-pd(3)+a34*pd(4)*cos(p(4));
xp43=pd(5)-b43*pd(6)*cos(p(6))-pd(3)+a43*pd(4)*cos(p(4));
xp35=pd(7)+c35*pd(8)*cos(p(8))-pd(3)+a35*pd(4)*cos(p(4));
xp53=pd(7)-c53*pd(8)*cos(p(8))-pd(3)+a53*pd(4)*cos(p(4));
xp56=pd(9)-pd(7)+c56*pd(8)*cos(p(8));
xp37=pd(10)-e37*pd(11)*cos(p(11))+ez2*pd(11)*sin(p(11))-pd(3)+a37*pd(4)*cos(p(4));


f13=(Feq13+P*A)*( (1+s*lnom13)/(1+s*(x13+lnom13)) )^k - P*A;
f23=(Feq23+P*A)*( (1+s*lnom23)/(1+s*(x23+lnom23)) )^k - P*A;



F10=k10*x10+d10*xp10-Feq10;
F13=f13+d13*xp13;
F34=k34*x34+d34*xp34-Feq34;
F35=k35*x35+d35*xp35-Feq35;
F56=k56*x56+d56*xp56-Feq56;
F20=k20*x20+d20*xp20-Feq20;
F23=f23+d23*xp23;
F43=k43*x43+d43*xp43-Feq43;
F53=k53*x53+d53*xp53-Feq53;
F37=k37*x37+d37*xp37-Feq37;

 
   M0=[m1 0  0  0  0  0  0  0  0  0  0
        0 m2 0  0  0  0  0  0  0  0  0
        0 0  m3 0  0  0  0  0  0  0  0
        0 0  0  l3 0  0  0  0  0  0  0
        0 0  0  0  m4 0  0  0  0  0  0
        0 0  0  0  0  l4 0  0  0  0  0
        0 0  0  0  0  0  m5 0  0  0  0
        0 0  0  0  0  0  0  l5 0  0  0
        0 0  0  0  0  0  0  0  m6 0  0
        0 0  0  0  0  0  0  0  0  m7 0
        0 0  0  0  0  0  0  0  0  0  l7];

    fa=[-F10+F13-m1*g
        -F20+F23-m2*g
        -F13-F23+F35+F34+F43+F53+F37-m3*g
        (a23*F23-a13*F13-a37*F37-a34*F34-a35*F35-a43*F43-a53*F53)*cos(p(4))
        -F43-F34-m4*g
        (b43*F43-b34*F34)*cos(p(6))
        -F53-F35+F56-m5*g
        (c53*F53-c35*F35-c56*F56)*cos(p(8))
        -F56-m6*g
        -F37-m7*g
        e37*cos(p(11))*F37];

    GT=[0
        0
        -1
        -az1*cos(p(4))+az2*sin(p(4))
        0
        0
        0
        0
        0
        1
        ez1*cos(p(11))+ez2*sin(p(11))];

    GTL = GT*lambda;

    g = -p(3)-az1*sin(p(4))-az2*cos(p(4))+p(10)+ez1*sin(p(11))-ez2*cos(p(11))+heq;
    
    f = [ M0*ddy(1:11) - fa + GTL
            g ];
        
  