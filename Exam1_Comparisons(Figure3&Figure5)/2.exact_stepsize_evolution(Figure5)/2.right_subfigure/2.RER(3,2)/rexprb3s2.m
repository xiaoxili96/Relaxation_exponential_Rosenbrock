function [exact_stepsize,tmesh]=rexprb3s2(tau)
options=optimset;  options = optimset(options,'Display','off');
T=4;  N=40;  h=2*pi/N;  xmesh=0:h:2*pi-h;  freq=-N/2:N/2-1;  area=2*pi;
[X,FREQ]=meshgrid(xmesh,freq);  W=(1/N)*exp(-1i*FREQ.*X);  W_hat=exp(1i*FREQ'.*X');
K=diag((1i*freq).^2);  d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  oc=zeros(2*d,1);  occ=zeros(d,1);
alpha=100;  func_F=@(x)alpha*(0.5*x.^2+0.25*x.^4);  func_f=@(x)alpha*(x+x.^3);  func_f_der=@(x)alpha*(1+3*x.^2);

tn=0;  Vn=zeros(d,1);  Un=1*(1+cos(xmesh))';  Vn_t=W*Vn;  Un_t=W*Un;  exact_stepsize=[];  tmesh=tn;
while (tn<(T-tau))
    Un=real(W_hat*Un_t);  energy_1=0.5*area*real(Vn_t'*Vn_t-Un_t'*K*Un_t);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);  
    pGpUn=[o -W*diag(func_f_der(Un))*W_hat;o o];  Jn=A+pGpUn;  G=[-W*func_f(Un);occ];  F=A*[Vn_t;Un_t]+G;
    %%%%  Un2
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
    exp_matrix=expm(Matrix);  Vm_t=Vn_t+exp_matrix(1:d,end);  Um_t=Un_t+exp_matrix(d+1:2*d,end);
    %%%%  Un+1
    Um=real(W_hat*Um_t);  Dn2=([-W*func_f(Um);occ]-G)-pGpUn*([Vm_t-Vn_t;Um_t-Un_t]);
    oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[Jn (2/tau/tau)*Dn2 oc F; oo II];
    exp_matrix=expm(Matrix);
    Vn_t_Update=exp_matrix(1:d,end);  Un_t_Update=exp_matrix(d+1:2*d,end);  Update_norm=norm(Vn_t_Update)+norm(Un_t_Update);
    if ( Update_norm==0 )
        gamma=1;
    else
        Un2=real(W_hat*Un_t_Update);  
        energy_2=area*real(Vn_t'*Vn_t_Update-Un_t'*K*Un_t_Update);
        energy_3=0.5*area*real(Vn_t_Update'*Vn_t_Update-Un_t_Update'*K*Un_t_Update);
        gamma=fsolve(@(gamma)func_of_gamma(Un,Un2,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma),1,options);
    end
    Vn_t_save=Vn_t;  Un_t_save=Un_t;  tn_save=tn;
    Vn_t=Vn_t+gamma*Vn_t_Update;  Un_t=Un_t+gamma*Un_t_Update;  tn=tn+gamma*tau
    tmesh=[tmesh tn];  exact_stepsize=[exact_stepsize gamma*tau];  
end

if ( (T-tn)<=0 )
    Vn_t=Vn_t_save;  Un_t=Un_t_save;  tn=tn_save;  tau=T-tn;  exact_stepsize=exact_stepsize(1:end-1);  tmesh=tmesh(1:end-1);
else
    tau=T-tn;
end
Un=real(W_hat*Un_t);  pGpUn=[o -W*diag(func_f_der(Un))*W_hat;o o];  Jn=A+pGpUn;  G=[-W*func_f(Un);occ];  F=A*[Vn_t;Un_t]+G;
%%%%  Un2
oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
exp_matrix=expm(Matrix);  Vm_t=Vn_t+exp_matrix(1:d,end);  Um_t=Un_t+exp_matrix(d+1:2*d,end);
%%%%  Un+1
Um=real(W_hat*Um_t);  Dn2=([-W*func_f(Um);occ]-G)-pGpUn*([Vm_t-Vn_t;Um_t-Un_t]);
oo=zeros(3,2*d);  II=[0 1 0;0 0 1;0 0 0];  Matrix=tau*[Jn (2/tau/tau)*Dn2 oc F; oo II];
exp_matrix=expm(Matrix); 
Vn_t=Vn_t+exp_matrix(1:d,end);  Un_t=Un_t+exp_matrix(d+1:2*d,end);  tn=tn+tau; 
