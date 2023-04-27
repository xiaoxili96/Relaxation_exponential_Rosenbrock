tic;
tau=1/200;
T=200;  N=40;  h=2*pi/N;  xmesh=0:h:2*pi-h;  freq=-N/2:N/2-1;  area=2*pi;
[X,FREQ]=meshgrid(xmesh,freq);  W=(1/N)*exp(-1i*FREQ.*X);  W_hat=exp(1i*FREQ'.*X');
K=diag((1i*freq).^2);  d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  oc=zeros(2*d,1);  occ=zeros(d,1);
alpha=100;  func_F=@(x)alpha*(0.5*x.^2+0.25*x.^4);  func_f=@(x)alpha*(x+x.^3);  func_f_der=@(x)alpha*(1+3*x.^2);

tn=0;  Vn=zeros(d,1);  Un=(1+cos(xmesh))';  Vn_t=W*Vn;  Un_t=W*Un;  GAMMA=[];  Energy=[];  tmesh=[];
while (tn<(T-0.5*tau))
    Un=real(W_hat*Un_t);  pGpUn=[o -W*diag(func_f_der(Un))*W_hat;o o];  
    Jn=A+pGpUn;  G=[-W*func_f(Un);occ];  F=A*[Vn_t;Un_t]+G;
    %%%%  Un2
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[(1/2)*Jn (1/2)*F; oo II];
    exp_matrix=expm(Matrix);  Vm_t=Vn_t+exp_matrix(1:d,end);  Um_t=Un_t+exp_matrix(d+1:2*d,end);
    %%%%  Un3
    Um=real(W_hat*Um_t);  Dn2=([-W*func_f(Um);occ]-G)-pGpUn*([Vm_t-Vn_t;Um_t-Un_t]);
    oo=zeros(1,2*d);  II=[0];  Matrix=tau*[Jn F; oo II];
    exp_matrix=expm(Matrix);  Vm_t=Vn_t+exp_matrix(1:d,end);  Um_t=Un_t+exp_matrix(d+1:2*d,end);
    %%%%  Un+1
    Um=real(W_hat*Um_t);  Dn3=([-W*func_f(Um);occ]-G)-pGpUn*([Vm_t-Vn_t;Um_t-Un_t]);
    oo=zeros(4,2*d);  II=[0 1 0 0;0 0 1 0;0 0 0 1;0 0 0 0];  
    Matrix=tau*[Jn (-48/tau/tau/tau)*Dn2+(12/tau/tau/tau)*Dn3 (16/tau/tau)*Dn2+(-2/tau/tau)*Dn3 oc F; oo II];
    energy_old=func_of_energy(Vn_t,Un_t,W_hat,K,area,h,func_F);  Energy=[Energy energy_old];  tmesh=[tmesh tn];
    exp_matrix=expm(Matrix);  Vn_t=Vn_t+exp_matrix(1:d,end);  Un_t=Un_t+exp_matrix(d+1:2*d,end);  tn=tn+tau
end
energy_old=func_of_energy(Vn_t,Un_t,W_hat,K,area,h,func_F);  Energy=[Energy energy_old];  tmesh=[tmesh tn];
toc;

Energy4s3=Energy;  tmesh4s3=tmesh;  time4s3=toc;
save('exprb4s3.mat','Energy4s3','tmesh4s3','time4s3');