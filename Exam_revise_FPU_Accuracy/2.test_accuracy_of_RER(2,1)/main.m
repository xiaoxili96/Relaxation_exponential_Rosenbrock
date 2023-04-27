load('reference.mat')
stepsize=[1/40 1/80 1/160 1/320];

yn_save=[];  xn_save=[];  Gamma_average_save=[];  err=[];
for k=1:size(stepsize,2)
    tau=stepsize(1,k);
    [yn,xn,gamma_average]=rexprb2s1(tau);
    yn_save=[yn_save yn];  
    xn_save=[xn_save xn];  
    Gamma_average_save=[Gamma_average_save gamma_average]; 
    err=[err max(abs([yn;xn]-Zn_c_100000))];
end


err 
Err_order=log(err(1:end-1)./err(2:end))./log(stepsize(1:end-1)./stepsize(2:end))
Gamma_average_save
Gamma_average_order=log(Gamma_average_save(1:end-1)./Gamma_average_save(2:end))./log(stepsize(1:end-1)./stepsize(2:end))