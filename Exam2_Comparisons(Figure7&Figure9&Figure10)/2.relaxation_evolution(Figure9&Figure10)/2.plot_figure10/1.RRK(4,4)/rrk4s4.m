function [GAMMA,tmesh]=rrk4s4(tau)
options=optimset;  options = optimset(options,'Display','off');
T=3;  N=20;  left=-7;  right=7;  bottom=-7;  top=7;  h=(right-left)/N;  
xmesh=left+0.5*h:h:right-0.5*h;  ymesh=xmesh;  [XMESH,YMESH]=meshgrid(xmesh,ymesh); 
KK=diag(ones(N-1,1),1)+diag(ones(N-1,1),-1);  I=eye(N);  K=kron(KK,I)+kron(I,KK);  K=K+(-1)*diag(sum(K,2));  K=(1/h/h)*K;
d=size(K,1);  o=zeros(d,d);  I=eye(d);  A=[o K;I o];  occ=zeros(d,1);
alpha=100;  func_F=@(x)alpha*(-cos(x));  func_f=@(x)alpha*(sin(x));  func_f_der=@(x)alpha*(cos(x));

tn=0;  Vn=zeros(d,1);  Un_temp=4*atan(exp(3-sqrt(XMESH.^2+YMESH.^2)));  Un=Un_temp(:);  GAMMA=[];  tmesh=tn;
while (tn<(T-tau))
    energy_1=0.5*h*h*real(Vn'*Vn-Un'*K*Un);  
    energy_old=func_of_energy(energy_1,Un,h,func_F);
    
    %%%% 1
    Vmid=Vn;  Umid=Un;
    Fn1=A*[Vmid;Umid]+[-func_f(Umid);occ];
    %%%% 2
    Vmid=Vn+(1/2)*tau*Fn1(1:N^2);  Umid=Un+(1/2)*tau*Fn1(N^2+1:2*N^2);
    Fn2=A*[Vmid;Umid]+[-func_f(Umid);occ];
    %%%% 3
    Vmid=Vn+((sqrt(2)-1)/2)*tau*Fn1(1:N^2)+((2-sqrt(2))/2)*tau*Fn2(1:N^2);  
    Umid=Un+((sqrt(2)-1)/2)*tau*Fn1(N^2+1:2*N^2)+((2-sqrt(2))/2)*tau*Fn2(N^2+1:2*N^2);
    Fn3=A*[Vmid;Umid]+[-func_f(Umid);occ];
    %%%% 4
    Vmid=Vn-(sqrt(2)/2)*tau*Fn2(1:N^2)+((2+sqrt(2))/2)*tau*Fn3(1:N^2);  
    Umid=Un-(sqrt(2)/2)*tau*Fn2(N^2+1:2*N^2)+((2+sqrt(2))/2)*tau*Fn3(N^2+1:2*N^2); 
    Fn4=A*[Vmid;Umid]+[-func_f(Umid);occ];
    %%%% next
    Vn1=Vn+(1/6)*tau*Fn1(1:N^2)+((2-sqrt(2))/6)*tau*Fn2(1:N^2)+((2+sqrt(2))/6)*tau*Fn3(1:N^2)+(1/6)*tau*Fn4(1:N^2);  
    Un1=Un+(1/6)*tau*Fn1(N^2+1:2*N^2)+((2-sqrt(2))/6)*tau*Fn2(N^2+1:2*N^2)+((2+sqrt(2))/6)*tau*Fn3(N^2+1:2*N^2)+(1/6)*tau*Fn4(N^2+1:2*N^2);
    Vn_Update=Vn1-Vn;  Un_Update=Un1-Un;
    Update_norm=norm(Vn_Update)+norm(Un_Update);
    if ( Update_norm==0 )
        gamma=1;
    else
        energy_2=h*h*real(Vn'*Vn_Update-Un'*K*Un_Update);
        energy_3=0.5*h*h*real(Vn_Update'*Vn_Update-Un_Update'*K*Un_Update);
        gamma=fsolve(@(gamma)func_of_gamma(Un,Un_Update,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma),1,options);
        fprintf('At time t=%d, the relaxation parameter is equal to %d\n',tn,gamma);
    end
    Vn_save=Vn;  Un_save=Un;  tn_save=tn;
    Vn=Vn+gamma*Vn_Update;  Un=Un+gamma*Un_Update;  tn=tn+gamma*tau;
    tmesh=[tmesh tn];  GAMMA=[GAMMA gamma];  
end

if ( (T-tn)<=0 )
    Vn=Vn_save;  Un=Un_save;  tn=tn_save;  tau=T-tn;  GAMMA=GAMMA(1:end-1);  tmesh=tmesh(1:end-1);
else
    tau=T-tn;
end
%%%% 1
Vmid=Vn;  Umid=Un;
Fn1=A*[Vmid;Umid]+[-func_f(Umid);occ];
%%%% 2
Vmid=Vn+(1/2)*tau*Fn1(1:N^2);  Umid=Un+(1/2)*tau*Fn1(N^2+1:2*N^2);
Fn2=A*[Vmid;Umid]+[-func_f(Umid);occ];
%%%% 3
Vmid=Vn+((sqrt(2)-1)/2)*tau*Fn1(1:N^2)+((2-sqrt(2))/2)*tau*Fn2(1:N^2);  
Umid=Un+((sqrt(2)-1)/2)*tau*Fn1(N^2+1:2*N^2)+((2-sqrt(2))/2)*tau*Fn2(N^2+1:2*N^2);
Fn3=A*[Vmid;Umid]+[-func_f(Umid);occ];
%%%% 4
Vmid=Vn-(sqrt(2)/2)*tau*Fn2(1:N^2)+((2+sqrt(2))/2)*tau*Fn3(1:N^2);  
Umid=Un-(sqrt(2)/2)*tau*Fn2(N^2+1:2*N^2)+((2+sqrt(2))/2)*tau*Fn3(N^2+1:2*N^2); 
Fn4=A*[Vmid;Umid]+[-func_f(Umid);occ];
%%%% next
Vn=Vn+(1/6)*tau*Fn1(1:N^2)+((2-sqrt(2))/6)*tau*Fn2(1:N^2)+((2+sqrt(2))/6)*tau*Fn3(1:N^2)+(1/6)*tau*Fn4(1:N^2);  
Un=Un+(1/6)*tau*Fn1(N^2+1:2*N^2)+((2-sqrt(2))/6)*tau*Fn2(N^2+1:2*N^2)+((2+sqrt(2))/6)*tau*Fn3(N^2+1:2*N^2)+(1/6)*tau*Fn4(N^2+1:2*N^2);
tn=tn+tau;