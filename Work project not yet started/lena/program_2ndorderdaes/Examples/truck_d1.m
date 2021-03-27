 function f = truck_d1(t,y,dy,ddy) 
 
 p=y(1:11);
lambda=10e-6*y(12);
pd=dy(1:11);
 
 
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


 
        
        f23_p3=(Feq23+P*A)*k*( (1+s*l)/(1+s*(x23+l)))^(k-1)*(-s*(1+s*l)/(1+s*(x23+l))^2);
        f23_p2=(Feq23+P*A)*k*( (1+s*l)/(1+s*(x23+l)))^(k-1)*(s*(1+s*l)/(1+s*(x23+l))^2);
        f23_p4=(Feq23+P*A)*k*( (1+s*l)/(1+s*(x23+l)))^(k-1)*(s*a23*cosp4*(1+s*l)/(1+s*(x23+l))^2);
        f13_p1=(Feq13+P*A)*k*( (1+s*l)/(1+s*(x13+l)))^(k-1)*(s*(1+s*l)/(1+s*(x13+l))^2);
        f13_p2=0;
        f13_p3=(Feq13+P*A)*k*( (1+s*l)/(1+s*(x13+l)))^(k-1)*(-s*(1+s*l)/(1+s*(x13+l))^2);
        f13_p4=(Feq13+P*A)*k*( (1+s*l)/(1+s*(x13+l)))^(k-1)*(-s*a13*cosp4*(1+s*l)/(1+s*(x13+l))^2);
        F13_p3=f13_p3;
        F13_p4=f13_p4-d13*a13*pd(4)*sinp4;
        F23_p2=f23_p2;
        F23_p4=f23_p4+d23*a23*pd(4)*sinp4;
        F37_p4= k37*a37*cosp4-d37*a37*pd(4)*sinp4;
        F37_p11=-k37*e37*cosp11+k37*ez2*sinp11+d37*e37*pd(11)*sinp11+d37*ez2*pd(11)*cosp11;
        F37_pd11=-d37*e37*cosp11+d37*ez2*sinp11;
        F37_pd4=d37*a37*cosp4;
        F34_p4= k34*a34*cosp4-d34*a34*pd(4)*sinp4;
        F34_p6=k34*b34*cosp6-d34*b34*pd(6)*sinp6;
        F34_pd4=d34*a34*cosp4;
        F34_pd6=d34*b34*cosp6;
        F35_p4=k35*a35*cosp4-d35*a35*pd(4)*sinp4;
        F35_p8=k35*c35*cosp8-d35*c35*pd(8)*sinp8;
        F35_pd4=d35*a35*cosp4;
        F35_pd8=d35*c35*cosp8;
        F53_p4=k53*a53*cosp4-d53*a53*pd(4)*sinp4;
        F53_p8=-k53*c53*cosp8+d53*c53*pd(8)*sinp8;
        F53_pd4=d53*a53*cosp4;
        F53_pd8=-d53*c53*cosp8;
        F43_p4= k43*a43*cosp4-d43*a43*pd(4)*sinp4;
        F43_p6=-k43*b43*cosp6+d43*b43*pd(6)*sinp6;
        F43_pd4=d43*a43*cosp4;
        F43_pd6=-d43*b43*cosp6;
        F56_p8=k56*c56*cosp8-d56*c56*pd(8)*sinp8;
        F56_pd8=d56*c56*cosp8;

        
        F1_p1  = k10-f13_p1;
        F1_p3  = -f13_p3;
        F1_p4  = -F13_p4;     
        F2_p2  = k20-f23_p2;
        F2_p3  = -f23_p3;
        F2_p4  = f23_p4-d23*a23*pd(4)*sinp4; 
        F3_p1  = f13_p1;
        F3_p2  = f23_p2;
        F3_p3  = f13_p3+f23_p3+k35+k34+k43+k53+k37;
        F3_p4  = f13_p4-d13*a13*pd(4)*sinp4+f23_p4+d23*a23*pd(4)*sinp4-k35*a35*cosp4+d35*a35*pd(4)*sinp4-k34*a34*cosp4+d34*a34*pd(4)*sinp4-k43*a43*cosp4+d43*a43*pd(4)*sinp4-k53*a53*cosp4+d53*a53*pd(4)*sinp4-k37*a37*cosp4+d37*a37*pd(4)*sinp4;
        F3_p5  = -k34-k43;
        F3_p6  = -k34*b34*cosp6+d34*b34*pd(6)*sinp6+k43*b43*cosp6-d43*b43*pd(6)*sinp6;
        F3_p7  = -k35-k53;
        F3_p8  = -k35*c35*cosp8+d35*c35*pd(8)*sinp8+k53*c53*cosp8-d53*c53*pd(8)*sinp8;    
        F3_p10 = -k37;
        F3_p11 = k37*e37*cosp11-k37*ez2*sinp11-d37*e37*pd(11)*sinp11-d37*ez2*pd(11)*cosp11;
        F4_p1  = a13*cosp4*f13_p1;
        F4_p2  = -a23*cosp4*f23_p2;
        F4_p3  = (-a23*f23_p3+a13*f13_p3-a37*k37-a34*k34-a35*k35-a43*k43-a53*k53)*cosp4;
        F4_p4  = (-a23*F23_p4+a13*F13_p4+a37*F37_p4+a34*F34_p4+a35*F35_p4+a43*F43_p4+a53*F53_p4)*cosp4-(-a23*F23+a13*F13+a37*F37+a34*F34+a34*F34+a35*F35+a43*F43+a53*F53)*sinp4+az1*lambda*sinp4+az2*cosp4*lambda;
        F4_p5  = (a34*k34+a43*k43)*cosp4;
        F4_p6  = (a34*F34_p6+a43*F43_p6)*cosp4;
        F4_p7  = (a35*k35+a53*k53)*cosp4;
        F4_p8  = (a35*F35_p8+a53*F53_p8)*cosp4; 
        F4_p10 = (a37*k37)*cosp4;
        F4_p11 = (a37*F37_p11)*cosp4;
        F4_L   = -az1*cosp4+az2*sinp4;   
        F5_p3  = -k43-k34;
        F5_p4  = k43*a43*cosp4-d43*a43*pd(4)*sinp4+k34*a34*cosp4-d34*a34*pd(4)*sinp4;
        F5_p5  = k43+k34;
        F5_p6  = -k43*b43*cosp6+d43*b43*pd(6)*sinp6+k34*b34*cosp6-d34*b34*pd(6)*sinp6;      
        F6_p3  = b43*cosp6*k43-b34*cosp6*k34;
        F6_p4  = -b43*cosp6*(k43*a43*cosp4-d43*a43*pd(4)*sinp4)+b34*cosp6*(k34*a34*cosp4-d34*a34*pd(4)*sinp4);
        F6_p5  = -b43*cosp6*k43+b34*cosp6*k34;
        F6_p6  = b43*sinp6*F43-b34*sinp6*F34-b43*cosp6*(-k43*b43*cosp6+d43*b43*pd(6)*sinp6)+b34*cosp6*(k34*b34*cosp6-d34*b34*pd(6)*sinp6);
        F7_p3  = -k53-k35;
        F7_p4  = k53*a53*cosp4-d53*a53*pd(4)*sinp4+k35*a35*cosp4-d35*a35*pd(4)*sinp4;    
        F7_p7  = k53+k35+k56;
        F7_p8  = -k53*c53*cosp8+d53*c53*pd(8)*sinp8+k35*c35*cosp8-d35*c35*pd(8)*sinp8-k56*c56*cosp8+d56*c56*pd(8)*sinp8;
        F7_p9  = -k56;
        F8_p3  = -cosp8*(-c53*k53+c35*k35);
        F8_p4  = -cosp8*(c53*F53_p4-c35*F35_p4);
        F8_p7  = -cosp8*(c53*k53-c35*k35+c56*k56);
        F8_p8  = sinp8*(c53*F53-c35*F35-c56*F56)-cosp8*(c53*F53_p8-c35*F35_p8-c56*F56_p8);
        F8_p9  = cosp8*(c56*k56);            
        F9_p7  = -k56;
        F9_p8  = k56*c56*cosp8-d56*c56*pd(8)*sinp8;
        F9_p9  = k56;
        F10_p3 = -k37;
        F10_p4 = k37*a37*cosp4-d37*a37*pd(4)*sinp4;    
        F10_p10= k37;
        F10_p11= -k37*e37*cosp11+k37*ez2*sinp11+d37*e37*pd(11)*sinp11+d37*ez2*pd(11)*cosp11;
        F10_L  = 1;            
        F11_p3 = k37*e37*cosp11;
        F11_p4 = -e37*cosp11*(k37*a37*cosp4-d37*a37*pd(4)*sinp4);     
        F11_p10= -e37*k37*cosp11;
        F11_p11= -e37*cosp11*(-k37*e37*cosp11+k37*ez2*sinp11+d37*e37*pd(11)*sinp11+d37*ez2*pd(11)*cosp11)+e37*sinp11*F37-ez1*lambda*sinp11+ez2*lambda*cosp11;
        F11_L  = ez1*cosp11+ez2*sinp11;
        g_p4   = -az1*cosp4+az2*sinp4;     
        g_p11  = ez1*cosp11+ez2*sinp11; 
                
            
        F1_pd1 = m1*a(1)+d10+d13;
        F1_pd3 = -d13;
        F1_pd4 = -d13*a13*cosp4;   
        F2_pd2 = m2*a(1)+d20+d23;
        F2_pd3 = -d23;
        F2_pd4 = d23*a23*cosp4;  
        F3_pd1 = -d13;
        F3_pd2 = -d23;
        F3_pd3 = m3*a(1)+d13+d23+d35+d34+d43+d53+d37;
        F3_pd4 = cosp4*(d13*a13-d23*a23-d35*a35-d34*a34-d43*a43-d53*a53-d37*a37);
        F3_pd5 = -d34-d43;
        F3_pd6 = -d34*b34*cosp6+d43*b43*cosp6;
        F3_pd7 = -d35-d53;
        F3_pd8 = -d35*c35*cosp8+d53*c53*cosp8;
        F3_pd10= -d37;
        F3_pd11= d37*e37*cosp11-d37*ez2*sinp11;
        F4_pd1 = cosp4*(-a13*d13);
        F4_pd2 = cosp4*(a23*d23);
        F4_pd3 = cosp4*(-a23*d23+a13*d13-a37*d37-a34*d34-a35*d35-a43*d43-a53*d53);
        F4_pd4 = l3*a(1)+cosp4*(a23*d23*a23*cosp4+a13*d13*a13*cosp4+a37*F37_pd4+a34*F34_pd4+a35*F35_pd4+a43*F43_pd4+a53*F53_pd4);
        F4_pd5 = cosp4*(a34*d34+a43*d43);
        F4_pd6 = cosp4*(a34*F34_pd6+a43*F43_pd6);
        F4_pd7 = cosp4*(a35*d35+a53*d53);
        F4_pd8 = cosp4*(a35*F35_pd8+a53*F53_pd8);
        F4_pd10= cosp4*(a37*d37);
        F4_pd11= cosp4*(a37*F37_pd11);
        F5_pd3 = -d43-d34;
        F5_pd4 = d43*a43*cosp4+d34*a34*cosp4;
        F5_pd5 = m4*a(1)+d43+d34;
        F5_pd6 = -d43*b43*cosp6+d34*b34*cosp6;    
        F6_pd3 = b43*cosp6*d43-b34*cosp6*d34;
        F6_pd4 = -b43*cosp6*d43*a43*cosp4+b34*cosp6*d34*a34*cosp4;
        F6_pd5 = -b43*cosp6*d43+b34*cosp6*d34;
        F6_pd6 = l4*a(1)+b43*cosp6*d43*b43*cosp6+b34*cosp6*d34*b34*cosp6;  
        F7_pd3 = -d53-d35;
        F7_pd4 = d53*a53*cosp4+d35*a35*cosp4;
        F7_pd7 = m5*a(1)+d53+d35+d56;
        F7_pd8 = -d53*c53*cosp8+d35*c35*cosp8-d56*c56*cosp8;
        F7_pd9 = -d56;    
        F8_pd3 = -cosp8*(-c53*d53+c35*d35);
        F8_pd4 = -cosp8*(c53*F53_pd4-c35*F35_pd4);
        F8_pd7 = -cosp8*(c53*d53-c35*k35+c56*d56);
        F8_pd8 = l5*a(1)-cosp8*(c53*F53_pd8-c35*F35_pd8-c56*F56_pd8);
        F8_pd9 = cosp8*c56*d56;    
        F9_pd7 = -d56;
        F9_pd8 = d56*c56*cosp8;
        F9_pd9 = m6*a(1)+d56;    
        F10_pd3= -d37;
        F10_pd4= d37*a37*cosp4;
        F10_pd10=m7*a(1)+d37;
        F10_pd11=-d37*e37*cosp11+d37*ez2*sinp11;
        F11_pd3= d37*e37*cosp11;
        F11_pd4= -d37*a37*cosp4*e37*cosp11;
        F11_pd10=-d37*e37*cosp11;
        F11_pd11=l7*a(1)+d37*e37*e37*cosp11^2-d37*ez2*sinp11*e37*cosp11;
        
        nq=11;
        

        

        M1=[F1_p1  0      F1_p3  F1_p4  0      0      0      0      0      0       0       0
            0      F2_p2  F2_p3  F2_p4  0      0      0      0      0      0       0       0
            F3_p1  F3_p2  F3_p3  F3_p4  F3_p5  F3_p6  F3_p7  F3_p8  0      F3_p10  F3_p11  -1
            F4_p1  F4_p2  F4_p3  F4_p4  F4_p5  F4_p6  F4_p7  F4_p8  0      F4_p10  F4_p11  F4_L
            0      0      F5_p3  F5_p4  F5_p5  F5_p6  0      0      0      0       0       0
            0      0      F6_p3  F6_p4  F6_p5  F6_p6  0      0      0      0       0       0
            0      0      F7_p3  F7_p4  0      0      F7_p7  F7_p8  F7_p9  0       0       0
            0      0      F8_p3  F8_p4  0      0      F8_p7  F8_p8  F8_p9  0       0       0
            0      0      0      0      0      0      F9_p7  F9_p8  F9_p9  0       0       0
            0      0      F10_p3 F10_p4 0      0      0      0      0      F10_p10 F10_p11 F10_L
            0      0      F11_p3 F11_p4 0      0      0      0      0      F11_p10 F11_p11 F11_L
            0      0     -1      g_p4   0      0      0      0      0      1       g_p11   0   ];
        
        
        M2=[F1_pd1  0       F1_pd3  F1_pd4  0       0       0       0       0       0        0       
            0       F2_pd2  F2_pd3  F2_pd4  0       0       0       0       0       0        0       
            F3_pd1  F3_pd2  F3_pd3  F3_pd4  F3_pd5  F3_pd6  F3_pd7  F3_pd8  0       F3_pd10  F3_pd11  
            F4_pd1  F4_pd2  F4_pd3  F4_pd4  F4_pd5  F4_pd6  F4_pd7  F4_pd8  0       F4_pd10  F4_pd11 
            0       0       F5_pd3  F5_pd4  F5_pd5  F5_pd6  0       0       0       0        0       
            0       0       F6_pd3  F6_pd4  F6_pd5  F6_pd6  0       0       0       0        0       
            0       0       F7_pd3  F7_pd4  0       0       F7_pd7  F7_pd8  F7_pd9  0        0       
            0       0       F8_pd3  F8_pd4  0       0       F8_pd7  F8_pd8  F8_pd9  0        0       
            0       0       0       0       0       0       F9_pd7  F9_pd8  F9_pd9  0        0       
            0       0       F10_pd3 F10_pd4 0       0       0       0       0       F10_pd10 F10_pd11 
            0       0       F11_pd3 F11_pd4 0       0       0       0       0       F11_pd10 F11_pd11 
            0       0       0       0       0       0       0       0       0       0        0      ];
        
            
        
        M=[ M1  M2];
        
        f = M1;
        