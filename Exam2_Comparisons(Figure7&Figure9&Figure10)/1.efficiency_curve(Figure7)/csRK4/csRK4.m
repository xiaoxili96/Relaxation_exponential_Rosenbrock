function [time,err,energy_err]=csRK4(tau)
tic;  
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
[coe1,coe2,coe3,coe4,coe5,coe6,GW]=generate_coefficient;

tn=0;  Vn=zeros(d,1);  Un_temp=4*atan(exp(3-sqrt(XMESH.^2+YMESH.^2)));  Un=Un_temp(:);  Zn=[Vn;Un];  ZZ=[Zn;Zn];  Energy=[];
while (tn<(T-0.5*tau))
    energy_1=0.5*h*h*real(Vn'*Vn-Un'*K*Un);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy_old];
    
    %%%%%%%%%
    Zn=ZZ(2*d+1:end,1);  Un=Zn(d+1:end,1);  Yn=A*Zn;  
    %%%%%%%%
    ZZ1=fsolve(@(ZZ1)equation(ZZ1,Zn,Un,Yn,func_f,A,d,tau,coe1,coe2,coe3,coe4,coe5,coe6,GW),ZZ,options);  ZZ=ZZ1;  Zn=ZZ1(2*d+1:end,1);  
    Vn_save=Zn(1:d,1);  Un_save=Zn(d+1:end,1);  tn_save=tn+tau;
    Vn=Vn_save;  Un=Un_save;  tn=tn_save
end
energy_1=0.5*h*h*real(Vn'*Vn-Un'*K*Un);  
energy_old=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy_old];
toc;

time=toc;
load('reference_1000_40.mat');  err=max(abs([Vn;Un]-[Vn_c_1000;Un_c_1000]));
energy_err=mean(abs(Energy-Energy(1))/abs(Energy(1)));