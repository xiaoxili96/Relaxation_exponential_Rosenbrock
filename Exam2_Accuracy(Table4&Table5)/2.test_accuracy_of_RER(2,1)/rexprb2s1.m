function [Vn,Un,gamma_average]=rexprb2s1(tau)
T=2;  N=40;  left=-7;  right=7;  bottom=-7; top=7;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 
KK=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  I=eye(N);  K=kron(KK,I)+kron(I,KK);  K=K+(-1)*diag(sum(K,2));  K=(1/h/h)*K;
d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  oc=zeros(2*d,1);  occ=zeros(d,1);
alpha=1;  func_F=@(x)alpha*(-cos(x));  func_f=@(x)alpha*(sin(x));  func_f_der=@(x)alpha*(cos(x));

tn=0;  Vn=zeros(d,1);  Un_temp=4*atan(exp(3-sqrt(XMESH.^2+YMESH.^2)));  Un=Un_temp(:);  GAMMA=[];
while (tn<(T-tau))
    pGpUn=[o -diag(func_f_der(Un));o o];  
    Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
    %%%%  Un+1
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
    exp_matrix=expm(Matrix);  Vn1=Vn+exp_matrix(1:d,end);  Un1=Un+exp_matrix(d+1:2*d,end);
    Vn_Update=Vn1-Vn;  Un_Update=Un1-Un;  
    energy_old=func_of_energy(Vn,Un,K,h,func_F);
    Update_norm=sum(abs(Vn_Update).^2)+sum(abs(Un_Update).^2);
    if ( Update_norm==0 )
        gamma=1;
    else
        gamma=fzero(@(gamma)func_of_gamma(Vn,Un,K,h,func_F,Vn_Update,Un_Update,energy_old,gamma),1);
        fprintf('distance=%d, Update=%d',abs(gamma-1),Update_norm);
    end
    Vn_save=Vn;  Un_save=Un;  tn_save=tn;
    GAMMA=[GAMMA gamma];  Vn=Vn+gamma*Vn_Update;  Un=Un+gamma*Un_Update;  tn=tn+gamma*tau
end

if ( (T-tn)<=0 )
    GAMMA=GAMMA(1:end-1);  
    Vn=Vn_save;  Un=Un_save;  tn=tn_save;  tau=T-tn;
else
    tau=T-tn;
end
pGpUn=[o -diag(func_f_der(Un));o o];  
Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
%%%%  Un+1
oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
exp_matrix=expm(Matrix);  Vn=Vn+exp_matrix(1:d,end);  Un=Un+exp_matrix(d+1:2*d,end);  tn=tn+tau
gamma_average=(sum(abs(GAMMA-1)))/(size(GAMMA,2));