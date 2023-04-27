function [Vn,Un]=exprb4s3(tau)
T=1;  N=10;  left=0;  right=1;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 
KK=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  I=eye(N);  K=kron(KK,I)+kron(I,KK);  K=K+(-1)*diag(sum(K,2));  K=(1/h/h)*K;
d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  oc=zeros(2*d,1);  occ=zeros(d,1);
alpha=100;  beta=10;  func_f=@(x)alpha*(sin(beta*x));  func_f_der=@(x)beta*alpha*(cos(beta*x));

tn=0;  Vn=zeros(d,1);  
Un_temp=sin(2*pi*XMESH).*sin(2*pi*YMESH);  
Un=Un_temp(:);  
while (tn<(T-0.5*tau))
    pGpUn=[o -diag(func_f_der(Un));o o];  
    Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
    %%%%  Un2
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[(1/2)*Jn (1/2)*F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%%  Un3
    Dn2=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%%  Un+1
    Dn3=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(4,2*d);  II=[0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];  
    Matrix=tau*[Jn (-48/tau/tau/tau)*Dn2+(12/tau/tau/tau)*Dn3 (16/tau/tau)*Dn2+(-2/tau/tau)*Dn3 oc F; oo II];
    exp_matrix=expm(Matrix);  Vn=Vn+exp_matrix(1:d,end);  Un=Un+exp_matrix(d+1:2*d,end);
    tn=tn+tau
end