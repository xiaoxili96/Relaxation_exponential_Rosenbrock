function [time,err,energy_err]=EAVF(tau)
tic;  warning('off');
options=optimset;
options = optimset(options,'TolX',1e-16);
options = optimset(options,'TolFun',1e-16);
options = optimset(options,'MaxFunEvals',Inf);
options = optimset(options,'MaxIter',2000);
options = optimset(options,'Display','off');
options = optimset(options,'Algorithm','levenberg-marquardt');
T=2;  N=40;  left=-7;  right=7;  bottom=-7; top=7;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 
KK=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  I=eye(N);  K=kron(KK,I)+kron(I,KK);  K=K+(-1)*diag(sum(K,2));  K=(1/h/h)*K;
d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o]; 
alpha=1;  func_F=@(x)alpha*(-cos(x));  func_f=@(x)alpha*(sin(x));
[GP,GW]=generate_GP_GW;  P0=expm(tau*A);  P1=phipade(tau*A,1);

tn=0;  Vn=zeros(d,1);  Un_temp=4*atan(exp(3-sqrt(XMESH.^2+YMESH.^2)));  Un=Un_temp(:);  Zn=[Vn;Un];  Energy=[];
while (tn<(T-0.5*tau))
    energy_1=0.5*h*h*real(Vn'*Vn-Un'*K*Un);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy_old];
    
    %%%%%%%%
    PZ=P0*Zn;
    %%%%%%%%
    Zn1=fsolve(@(Zn1)equation(Zn1,Un,PZ,P1,d,tau,func_f,GP,GW),Zn,options);  Zn=Zn1;  
    Vn_save=Zn(1:d,1);  Un_save=Zn(d+1:end,1);  tn_save=tn+tau;
    Vn=Vn_save;  Un=Un_save;  tn=tn_save
end
energy_1=0.5*h*h*real(Vn'*Vn-Un'*K*Un);  
energy_old=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy_old];
toc;

time=toc;
load('reference_1000_40.mat');  err=max(abs([Vn;Un]-[Vn_c_1000;Un_c_1000]));
energy_err=mean(abs(Energy-Energy(1))/abs(Energy(1)));