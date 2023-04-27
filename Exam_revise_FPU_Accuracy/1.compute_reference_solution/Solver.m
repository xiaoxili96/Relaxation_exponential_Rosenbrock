function [Zn,t]=Solver(tau)
T=1;
%%%%%%%%%%%%%%%% 2-stages %%%%%%%%%%%%%%%%
% Gauss
% A=[1/4 1/4-sqrt(3)/6; 1/4+sqrt(3)/6 1/4];
% b=[1/2 1/2];
% c=[1/2-sqrt(3)/6; 1/2+sqrt(3)/6];
%%%%%%%%%%%%%%%% 3-stages %%%%%%%%%%%%%%%%
% Gauss
% A=[5/36 2/9-sqrt(15)/15 5/36-sqrt(15)/30; ...
%    5/36+sqrt(15)/24 2/9 5/36-sqrt(15)/24; ...
%    5/36+sqrt(15)/30 2/9+sqrt(15)/15 5/36];
% b=[5/18 4/9 5/18];
% c=[1/2-sqrt(15)/10; 1/2; 1/2+sqrt(15)/10];
%%%%%%%%%%%%%%%% 4-stages %%%%%%%%%%%%%%%%
% Gauss
w1=1/8-sqrt(30)/144;  w2=0.5*sqrt((15+2*sqrt(30))/35);  w3=w2*(1/6+sqrt(30)/24); 
w4=w2*(1/21+5*sqrt(30)/168);  w5=w2-2*w3;
w11=1/8+sqrt(30)/144;  w22=0.5*sqrt((15-2*sqrt(30))/35);  w33=w22*(1/6-sqrt(30)/24); 
w44=w22*(1/21-5*sqrt(30)/168);  w55=w22-2*w33;
A=[w1 w11-w3+w44 w11-w3-w44 w1-w5; ... 
   w1-w33+w4 w11 w11-w55 w1-w33-w4; ...
   w1+w33+w4 w11+w55 w11 w1+w33-w4; ...
   w1+w5 w11+w3+w44 w11+w3-w44 w1];
b=[2*w1 2*w11 2*w11 2*w1];
c=[0.5-w2; 0.5-w22; 0.5+w22; 0.5+w2];

m=3;  w=50;  om=zeros(m,m);  Im=eye(m);  o2m=zeros(2*m,2*m);  I2m=eye(2*m);
AA=[om om; om w^2*Im];  L=[o2m -AA; I2m o2m];
d=4*m;  s=size(b,2);  Id=eye(d);  es=ones(s,1);  Ids=eye(s*d,s*d); 
k_A_Id=kron(A,Id);  k_A_L=kron(A,L);  k_b_Id=kron(b,Id);  k_b_L=kron(b,L);
Zn=[1;0;0;1;0;0;1;0;0;1/w;0;0];  t=0; 

for k=1:round(T/tau)
    k_es_Zn=kron(es,Zn);  Zmid=k_es_Zn;  Iter_err=1;  count=0; 
    while (Iter_err>10^(-16) && count < 100)
        F=compute_nonlinear(Zmid,d,s,m);
        Vector=Zmid-k_es_Zn-tau*k_A_L*Zmid-tau*k_A_Id*F; 
        F_der=compute_nonlinear_der(Zmid,d,s,m);
        Matrix=Ids-tau*k_A_L-tau*k_A_Id*F_der;
        Zmid_save=Zmid;
        Zmid=Zmid-Matrix\Vector;
        Iter_err=max(abs(Zmid-Zmid_save));
        count=count+1;
    end
    F=compute_nonlinear(Zmid,d,s,m);
    Zn=Zn+tau*k_b_L*Zmid+tau*k_b_Id*F;
    t=t+tau; 
    k
end