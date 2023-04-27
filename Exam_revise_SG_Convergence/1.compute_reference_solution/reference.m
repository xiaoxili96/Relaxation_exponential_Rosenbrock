tau=1/100000;
T=1;  N=10;  left=0;  right=1;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 
KK=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  I=eye(N);  K=kron(KK,I)+kron(I,KK);  K=K+(-1)*diag(sum(K,2));  K=(1/h/h)*K;
d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  oc=zeros(2*d,1);  occ=zeros(d,1);
alpha=100;  beta=10;  func_F=@(x)(1/beta)*alpha*(-cos(beta*x));  func_f=@(x)alpha*(sin(beta*x));  func_f_der=@(x)beta*alpha*(cos(beta*x));

tn=0;  Vn=zeros(d,1);  
Un_temp=sin(2*pi*XMESH).*sin(2*pi*YMESH); 
Un=Un_temp(:);
while (tn<(T-tau))
    pGpUn=[o -diag(func_f_der(Un));o o];  
    Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
    %%%%  Un2
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[(1/2)*Jn (1/2)*F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%%  Un3
    Dn2=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[(1/2)*Jn (1/tau/tau)*Dn2 oc (1/2)*F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%% Un4
    Dn3=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[(1/3)*Jn (8/27/tau/tau)*Dn3 oc (1/3)*F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%% Un5
    Dn4=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[Jn (18/tau/tau)*Dn4 oc F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%% Un+1
    Dn5=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(5,2*d);  II=[0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 0];   
    Matrix=tau*[Jn (-1152/tau/tau/tau/tau)*Dn3+(1944/tau/tau/tau/tau)*Dn4+(72/tau/tau/tau/tau)*Dn5 (384/tau/tau/tau)*Dn3+(-729/tau/tau/tau)*Dn4+(-15/tau/tau/tau)*Dn5 (-32/tau/tau)*Dn3+(81/tau/tau)*Dn4+(1/tau/tau)*Dn5 oc F; oo II];
    exp_matrix=expm(Matrix);  Vn1=Vn+exp_matrix(1:d,end);  Un1=Un+exp_matrix(d+1:2*d,end);
    Vn_Update=Vn1-Vn;  Un_Update=Un1-Un;  
    energy_old=func_of_energy(Vn,Un,K,h,func_F);
    Update_norm=sum(abs(Vn_Update).^2)+sum(abs(Un_Update).^2);
    if ( Update_norm==0 )
        gamma=1;
    else
        gamma=fzero(@(gamma)func_of_gamma(Vn,Un,K,h,func_F,Vn_Update,Un_Update,energy_old,gamma),1); 
    end
    Vn_save=Vn;  Un_save=Un;  tn_save=tn;
    Vn=Vn+gamma*Vn_Update;  Un=Un+gamma*Un_Update;  tn=tn+gamma*tau
end

if ( (T-tn)<=0 )
    Vn=Vn_save;  Un=Un_save;  tn=tn_save;  tau=T-tn;
else
    tau=T-tn;
end
pGpUn=[o -diag(func_f_der(Un));o o];  
Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
%%%%  Un2
oo=zeros(1,2*d);  II=[0];  Matrix=tau*[(1/2)*Jn (1/2)*F; oo II];
exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
%%%%  Un3
Dn2=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[(1/2)*Jn (1/tau/tau)*Dn2 oc (1/2)*F; oo II];
exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
%%%% Un4
Dn3=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[(1/3)*Jn (8/27/tau/tau)*Dn3 oc (1/3)*F; oo II];
exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
%%%% Un5
Dn4=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[Jn (18/tau/tau)*Dn4 oc F; oo II];
exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
%%%% Un+1
Dn5=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
oo=zeros(5,2*d);  II=[0 1 0 0 0;0 0 1 0 0;0 0 0 1 0;0 0 0 0 1;0 0 0 0 0];   
Matrix=tau*[Jn (-1152/tau/tau/tau/tau)*Dn3+(1944/tau/tau/tau/tau)*Dn4+(72/tau/tau/tau/tau)*Dn5 (384/tau/tau/tau)*Dn3+(-729/tau/tau/tau)*Dn4+(-15/tau/tau/tau)*Dn5 (-32/tau/tau)*Dn3+(81/tau/tau)*Dn4+(1/tau/tau)*Dn5 oc F; oo II];
exp_matrix=expm(Matrix);  Vn=Vn+exp_matrix(1:d,end);  Un=Un+exp_matrix(d+1:2*d,end);  tn=tn+tau

Vn_c_100000=Vn;  Un_c_100000=Un;  save('reference.mat','Vn_c_100000','Un_c_100000');