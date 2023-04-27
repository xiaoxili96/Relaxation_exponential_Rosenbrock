function [exact_stepsize,tmesh]=rrk3s3(tau)
options=optimset;  options = optimset(options,'Display','off');
T=3;  N=40;  h=2*pi/N;  xmesh=0:h:2*pi-h;  freq=-N/2:N/2-1;  area=2*pi;
[X,FREQ]=meshgrid(xmesh,freq);  W=(1/N)*exp(-1i*FREQ.*X);  W_hat=exp(1i*FREQ'.*X');
K=diag((1i*freq).^2);  d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  occ=zeros(d,1);
alpha=100;  func_F=@(x)alpha*(0.5*x.^2+0.25*x.^4);  func_f=@(x)alpha*(x+x.^3);

tn=0;  Vn=zeros(d,1);  Un=1*(1+cos(xmesh))';  Vn_t=W*Vn;  Un_t=W*Un;  exact_stepsize=[];  tmesh=tn;
while (tn<(T-tau))
    Un=real(W_hat*Un_t);  energy_1=0.5*area*real(Vn_t'*Vn_t-Un_t'*K*Un_t);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);
    %%%% 1
    Vmid_t=Vn_t;  Umid_t=Un_t;
    Umid=real(W_hat*Umid_t);  Fn1=A*[Vmid_t;Umid_t]+[-W*func_f(Umid);occ];
    %%%% 2
    Vmid_t=Vn_t+(1/3)*tau*Fn1(1:N);  Umid_t=Un_t+(1/3)*tau*Fn1(N+1:2*N);
    Umid=real(W_hat*Umid_t);  Fn2=A*[Vmid_t;Umid_t]+[-W*func_f(Umid);occ];
    %%%% 3
    Vmid_t=Vn_t+(2/3)*tau*Fn2(1:N);  Umid_t=Un_t+(2/3)*tau*Fn2(N+1:2*N);
    Umid=real(W_hat*Umid_t);  Fn3=A*[Vmid_t;Umid_t]+[-W*func_f(Umid);occ];
    %%%% next
    Vn1_t=Vn_t+(1/4)*tau*Fn1(1:N)+(3/4)*tau*Fn3(1:N);  Un1_t=Un_t+(1/4)*tau*Fn1(N+1:2*N)+(3/4)*tau*Fn3(N+1:2*N);
    Vn_t_Update=Vn1_t-Vn_t;  Un_t_Update=Un1_t-Un_t;
    Update_norm=norm(Vn_t_Update)+norm(Un_t_Update);
    if ( Update_norm==0 )
        gamma=1;
    else
        Un2=real(W_hat*Un_t_Update);  
        energy_2=area*real(Vn_t'*Vn_t_Update-Un_t'*K*Un_t_Update);
        energy_3=0.5*area*real(Vn_t_Update'*Vn_t_Update-Un_t_Update'*K*Un_t_Update);
        gamma=fsolve(@(gamma)func_of_gamma(Un,Un2,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma),1,options);
        fprintf('At time t=%d, the exact stepsize is equal to %d\n',tn,gamma*tau);
    end
    Vn_t_save=Vn_t;  Un_t_save=Un_t;  tn_save=tn;
    Vn_t=Vn_t+gamma*Vn_t_Update;  Un_t=Un_t+gamma*Un_t_Update;  tn=tn+gamma*tau;
    tmesh=[tmesh tn];  exact_stepsize=[exact_stepsize gamma*tau];  
end

if ( (T-tn)<=0 )
    Vn_t=Vn_t_save;  Un_t=Un_t_save;  tn=tn_save;  tau=T-tn;  exact_stepsize=exact_stepsize(1:end-1);  tmesh=tmesh(1:end-1);
else
    tau=T-tn;
end
%%%% 1
Vmid_t=Vn_t;  Umid_t=Un_t;
Umid=real(W_hat*Umid_t);  Fn1=A*[Vmid_t;Umid_t]+[-W*func_f(Umid);occ];
%%%% 2
Vmid_t=Vn_t+(1/3)*tau*Fn1(1:N);  Umid_t=Un_t+(1/3)*tau*Fn1(N+1:2*N);
Umid=real(W_hat*Umid_t);  Fn2=A*[Vmid_t;Umid_t]+[-W*func_f(Umid);occ];
%%%% 3
Vmid_t=Vn_t+(2/3)*tau*Fn2(1:N);  Umid_t=Un_t+(2/3)*tau*Fn2(N+1:2*N);
Umid=real(W_hat*Umid_t);  Fn3=A*[Vmid_t;Umid_t]+[-W*func_f(Umid);occ];
%%%% next
Vn_t=Vn_t+(1/4)*tau*Fn1(1:N)+(3/4)*tau*Fn3(1:N);  Un_t=Un_t+(1/4)*tau*Fn1(N+1:2*N)+(3/4)*tau*Fn3(N+1:2*N);  
tn=tn+tau;
