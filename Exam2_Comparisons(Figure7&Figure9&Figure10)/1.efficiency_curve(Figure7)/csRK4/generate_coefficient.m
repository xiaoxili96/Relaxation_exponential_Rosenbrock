function [coe1,coe2,coe3,coe4,coe5,coe6,GW]=generate_coefficient

p1=-0.9061798459386640;  p2=-0.5384693101056830;  p3=0;  p4=-p2;  p5=-p1;  GP=[p1 p2 p3 p4 p5];

w1=0.2369268850561890;  w2=0.4786286704993660;  w3=0.5688888888888880;  w4=w2;  w5=w1;  GW=[w1 w2 w3 w4 w5];

GP=(1/2)*(GP+1);  GW=(1/2)*GW;

c1=0;  c2=1/2;  c3=1;
coe1=2*(GP-c2).*(GP-c3);
coe2=(-4)*GP.*(GP-c3);
coe3=2*(GP-c1).*(GP-c2);

coe4=((5/4)-(3/2)*GP).*GW;

syms x;
l=((5/4)-(3/2)*x);
a=2*(x-(1/2))*(x-1);
b=(-4)*x*(x-1);
c=2*x*(x-(1/2));

a1=double(int(a*l,0,1));
a2=double(int(b*l,0,1));
a3=double(int(c*l,0,1));  coe5=[a1 a2 a3];
aa1=double(int(a,0,1));
aa2=double(int(b,0,1));
aa3=double(int(c,0,1));  coe6=[aa1 aa2 aa3];