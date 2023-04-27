function [time,err,energy_err]=rexprb2s1(tau)
tic;  
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

tn=0;  Vn=zeros(d,1);  Un=1*(1+cos(xmesh))';  Vn_t=W*Vn;  Un_t=W*Un;  Energy=[];
while (tn<(T-tau))
    Un=real(W_hat*Un_t);  energy_1=0.5*area*real(Vn_t'*Vn_t-Un_t'*K*Un_t);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);  Energy=[Energy energy_old];
    
    pGpUn=[o -W*diag(func_f_der(Un))*W_hat;o o];  Jn=A+pGpUn;  G=[-W*func_f(Un);occ];  F=A*[Vn_t;Un_t]+G;
    %%%%  Un+1
    oo=zeros(1,2*d);  II=0;  Matrix=tau*[Jn F; oo II];  exp_matrix=expm(Matrix);
    Vn_t_Update=exp_matrix(1:d,end);  Un_t_Update=exp_matrix(d+1:2*d,end);  Update_norm=norm(Vn_t_Update)+norm(Un_t_Update);
    if ( Update_norm==0 )
        gamma=1;
    else
        Un2=real(W_hat*Un_t_Update);  
        energy_2=area*real(Vn_t'*Vn_t_Update-Un_t'*K*Un_t_Update);
        energy_3=0.5*area*real(Vn_t_Update'*Vn_t_Update-Un_t_Update'*K*Un_t_Update);
        gamma=fsolve(@(gamma)func_of_gamma(Un,Un2,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma),1,options);
        % fprintf('distance=%d, Update=%d\n',abs(gamma-1),Update_norm);
    end
    Vn_t_save=Vn_t;  Un_t_save=Un_t;  tn_save=tn;
    Vn_t=Vn_t+gamma*Vn_t_Update;  Un_t=Un_t+gamma*Un_t_Update;  tn=tn+gamma*tau;
end

if ( (T-tn)<=0 )
    Vn_t=Vn_t_save;  Un_t=Un_t_save;  tn=tn_save;  tau=T-tn;
else
    tau=T-tn;
end
Un=real(W_hat*Un_t);  pGpUn=[o -W*diag(func_f_der(Un))*W_hat;o o];  
Jn=A+pGpUn;  G=[-W*func_f(Un);occ];  F=A*[Vn_t;Un_t]+G;
%%%%  Un+1
oo=zeros(1,2*d);  II=0;  Matrix=tau*[Jn F; oo II];  exp_matrix=expm(Matrix);  
Vn_t=Vn_t+exp_matrix(1:d,end);  Un_t=Un_t+exp_matrix(d+1:2*d,end);  tn=tn+tau;
toc;

time=toc;
load('reference.mat');  err=max(abs([Vn_t;Un_t]-[Vn_t_c;Un_t_c]));
energy_err=mean(abs(Energy-Energy(1))/abs(Energy(1)));