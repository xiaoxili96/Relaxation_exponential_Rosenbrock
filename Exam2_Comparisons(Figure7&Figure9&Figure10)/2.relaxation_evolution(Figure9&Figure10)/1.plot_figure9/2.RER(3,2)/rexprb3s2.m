tau=1/16;
options=optimset;  options = optimset(options,'Display','off');
T=3;  N=20;  left=-7;  right=7;  bottom=-7;  top=7;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 
KK=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  I=eye(N);  K=kron(KK,I)+kron(I,KK);  K=K+(-1)*diag(sum(K,2));  K=(1/h/h)*K;
d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  oc=zeros(2*d,1);  occ=zeros(d,1);
alpha=1000;  func_F=@(x)alpha*(-cos(x));  func_f=@(x)alpha*(sin(x));  func_f_der=@(x)alpha*(cos(x));

tn=0;  Vn=zeros(d,1);  Un_temp=4*atan(exp(3-sqrt(XMESH.^2+YMESH.^2))); 
Un=Un_temp(:);  GAMMA=[];  tmesh=tn;  Energy=[];
while (tn<(T-tau))
    energy_1=0.5*h*h*real(Vn'*Vn-Un'*K*Un);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy_old];
    
    pGpUn=[o -diag(func_f_der(Un));o o];  
    Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
    %%%%  Un2
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
    exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
    %%%%  Un+1
    Dn2=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
    oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[Jn (2/tau/tau)*Dn2 oc F; oo II];  exp_matrix=expm(Matrix);
    Vn_Update=exp_matrix(1:d,end);  Un_Update=exp_matrix(d+1:2*d,end);  Update_norm=norm(Vn_Update)+norm(Un_Update);
    if ( Update_norm==0 )
        gamma=1;
    else
        energy_2=h*h*real(Vn'*Vn_Update-Un'*K*Un_Update);
        energy_3=0.5*h*h*real(Vn_Update'*Vn_Update-Un_Update'*K*Un_Update);
        gamma=fsolve(@(gamma)func_of_gamma(Un,Un_Update,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma),1,options);
    end
    Vn_save=Vn;  Un_save=Un;  tn_save=tn;
    Vn=Vn+gamma*Vn_Update;  Un=Un+gamma*Un_Update;  tn=tn+gamma*tau
    tmesh=[tmesh tn];  GAMMA=[GAMMA gamma];  
end

if ( (T-tn)<=0 )
    Vn=Vn_save;  Un=Un_save;  tn=tn_save;  tau=T-tn;  GAMMA=GAMMA(1:end-1);  tmesh=tmesh(1:end-1);
else
    tau=T-tn;
end
pGpUn=[o -diag(func_f_der(Un));o o];  
Jn=A+pGpUn;  G=[-func_f(Un);occ];  F=A*[Vn;Un]+G;
%%%%  Un2
oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
exp_matrix=expm(Matrix);  Vm=Vn+exp_matrix(1:d,end);  Um=Un+exp_matrix(d+1:2*d,end);
%%%%  Un+1
Dn2=([-func_f(Um);occ]-G)-pGpUn*([Vm-Vn;Um-Un]);
oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[Jn (2/tau/tau)*Dn2 oc F; oo II];  exp_matrix=expm(Matrix); 
Vn=Vn+exp_matrix(1:d,end);  Un=Un+exp_matrix(d+1:2*d,end);  tn=tn+tau;

save('data_rexprb3s2.mat','GAMMA','tmesh');