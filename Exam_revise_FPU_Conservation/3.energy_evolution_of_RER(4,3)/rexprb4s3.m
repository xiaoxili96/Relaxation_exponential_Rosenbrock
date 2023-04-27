tau=1/4;
T=300;  m=3;  w=50;  om=zeros(m,m);  Im=eye(m);  o2m=zeros(2*m,2*m);  I2m=eye(2*m);  
AA=[om om; om w^2*Im];  A=[o2m -AA; I2m o2m];  d=4*m;  oc=zeros(d,1); 
Zn=[1;0;0;1;0;0;1;0;0;1/w;0;0];
yn=Zn(1:2*m);  xn=Zn(2*m+1:4*m);  tn=0;
Energy=[];  tmesh=[];  GAMMA=[];

while (tn<(T-tau))
    pGpUn=compute_nonlinear_der(xn,m); 
    Jn=A+pGpUn;  G=compute_nonlinear(xn,m);  F=A*[yn;xn]+G;
    %%%%  Un2
    oo=zeros(1,d);  II=[0];  Matrix=tau*[(1/2)*Jn (1/2)*F; oo II];
    exp_matrix=expm(Matrix);  ym=yn+exp_matrix(1:2*m,end);  xm=xn+exp_matrix(2*m+1:d,end);
    %%%%  Un3
    Dn2=(compute_nonlinear(xm,m)-G)-pGpUn*([ym-yn;xm-xn]);
    oo=zeros(1,d);  II=[0];  Matrix=tau*[Jn F; oo II];
    exp_matrix=expm(Matrix);  ym=yn+exp_matrix(1:2*m,end);  xm=xn+exp_matrix(2*m+1:d,end);
    %%%%  Un+1
    Dn3=(compute_nonlinear(xm,m)-G)-pGpUn*([ym-yn;xm-xn]);
    oo=zeros(4,d);  II=[0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];  
    Matrix=tau*[Jn (-48/tau/tau/tau)*Dn2+(12/tau/tau/tau)*Dn3 (16/tau/tau)*Dn2+(-2/tau/tau)*Dn3 oc F; oo II];
    exp_matrix=expm(Matrix);  yn1=yn+exp_matrix(1:2*m,end);  xn1=xn+exp_matrix(2*m+1:d,end);
    yn_Update=yn1-yn;  xn_Update=xn1-xn;  
    energy_old=func_of_energy(yn,xn,AA,m);
    Energy=[Energy energy_old];  tmesh=[tmesh tn]; 
    Update_norm=sum(abs(yn_Update).^2)+sum(abs(xn_Update).^2);
    if ( Update_norm==0 )
        gamma=1;
    else
        gamma=fzero(@(gamma)func_of_gamma(yn,xn,AA,m,yn_Update,xn_Update,energy_old,gamma),1);
    end
    yn_save=yn;  xn_save=xn;  tn_save=tn;
    GAMMA=[GAMMA gamma];  yn=yn+gamma*yn_Update;  xn=xn+gamma*xn_Update;  tn=tn+gamma*tau
end
energy_old=func_of_energy(yn,xn,AA,m);  Energy=[Energy energy_old];  tmesh=[tmesh tn];

if ( (T-tn)<=0 )
    GAMMA=GAMMA(1:end-1);  Energy=Energy(1:end-1);  tmesh=tmesh(1:end-1);
    yn=yn_save;  xn=xn_save;  tn=tn_save;  tau=T-tn;
else
        tau=T-tn;
end
pGpUn=compute_nonlinear_der(xn,m); 
Jn=A+pGpUn;  G=compute_nonlinear(xn,m);  F=A*[yn;xn]+G;
%%%%  Un2
oo=zeros(1,d);  II=[0];  Matrix=tau*[(1/2)*Jn (1/2)*F; oo II];
exp_matrix=expm(Matrix);  ym=yn+exp_matrix(1:2*m,end);  xm=xn+exp_matrix(2*m+1:d,end);
%%%%  Un3
Dn2=(compute_nonlinear(xm,m)-G)-pGpUn*([ym-yn;xm-xn]);
oo=zeros(1,d);  II=[0];  Matrix=tau*[Jn F; oo II];
exp_matrix=expm(Matrix);  ym=yn+exp_matrix(1:2*m,end);  xm=xn+exp_matrix(2*m+1:d,end);
%%%%  Un+1
Dn3=(compute_nonlinear(xm,m)-G)-pGpUn*([ym-yn;xm-xn]);
oo=zeros(4,d);  II=[0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];  
Matrix=tau*[Jn (-48/tau/tau/tau)*Dn2+(12/tau/tau/tau)*Dn3 (16/tau/tau)*Dn2+(-2/tau/tau)*Dn3 oc F; oo II];
exp_matrix=expm(Matrix);  yn=yn+exp_matrix(1:2*m,end);  xn=xn+exp_matrix(2*m+1:d,end);  tn=tn+tau
gamma_average=(sum(abs(GAMMA-1)))/(size(GAMMA,2));

plot(tmesh,abs(Energy-Energy(1))/abs(Energy(1)));
Energyr4s3=Energy;  tmeshr4s3=tmesh;
save('rexprb4s3.mat','Energyr4s3','tmeshr4s3');