function [time,err,energy_err]=EAVF(tau)
tic;  warning('off');
options=optimset;
options = optimset(options,'TolX',1e-16);
options = optimset(options,'TolFun',1e-16);
options = optimset(options,'MaxFunEvals',Inf);
options = optimset(options,'MaxIter',2000);
options = optimset(options,'Display','off');
options = optimset(options,'Algorithm','levenberg-marquardt');
T=1;  N=40;  h=2*pi/N;  xmesh=0:h:2*pi-h;  freq=-N/2:N/2-1;  area=2*pi;
[X,FREQ]=meshgrid(xmesh,freq);  W=(1/N)*exp(-1i*FREQ.*X);  W_hat=exp(1i*FREQ'.*X');
K=diag((1i*freq).^2);  d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  occ=zeros(d,1);
alpha=100;  func_F=@(x)alpha*(0.5*x.^2+0.25*x.^4);  func_f=@(x)alpha*(x+x.^3);  func_f_der=@(x)alpha*(1+3*x.^2);

tn=0;  Vn=zeros(d,1);  Un=1*(1+cos(xmesh))';  Vn_t=W*Vn;  Un_t=W*Un;  Zn_t=[Vn_t;Un_t];  Energy=[];
P0=expm(tau*A);  P1=phipade(tau*A,1);
while (tn<(T-0.5*tau))
    Un=real(W_hat*Un_t);  energy_1=0.5*area*real(Vn_t'*Vn_t-Un_t'*K*Un_t);  
    energy=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy];
    
    %%%%%%%%
    Un_t=Zn_t(d+1:end,1);  old=real(W_hat*Un_t);  old1=old+old.^3;  PZ=P0*Zn_t;
    %%%%%%%%
    Zn1_t=fsolve(@(Zn1_t)equation(Zn1_t,W,W_hat,PZ,P1,d,tau,alpha,old,old1),Zn_t,options);  Zn_t=Zn1_t;  
    Vn_t_save=Zn_t(1:d,1);  Un_t_save=Zn_t(d+1:end,1);  tn_save=tn+tau;
    Vn_t=Vn_t_save;  Un_t=Un_t_save;  tn=tn_save;
end
Un=real(W_hat*Un_t);  energy_1=0.5*area*real(Vn_t'*Vn_t-Un_t'*K*Un_t);  
energy=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy];
toc;

time=toc;
load('reference.mat');  err=max(abs([Vn_t;Un_t]-[Vn_t_c;Un_t_c]));
energy_err=mean(abs(Energy-Energy(1))/abs(Energy(1)));