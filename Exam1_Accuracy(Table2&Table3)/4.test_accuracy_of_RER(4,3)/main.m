load('reference.mat')
stepsize=[1/600 1/1200 1/1800 1/2400];

V_save=[];  U_save=[];  Gamma_average_save=[];
for k=1:size(stepsize,2)
    tau=stepsize(1,k);
    [Vn_t,Un_t,gamma_average]=rexprb4s3(tau);
    V_save=[V_save Vn_t];  
    U_save=[U_save Un_t];  
    Gamma_average_save=[Gamma_average_save gamma_average]; 
end

err=[];
for k=1:(size(U_save,2))
    err=[err max([abs(V_save(:,k)-Vn_t_c_10000);abs(U_save(:,k)-Un_t_c_10000)])];
end

err
Err_order=log(err(1:end-1)./err(2:end))./log(stepsize(1:end-1)./stepsize(2:end))
Gamma_average_save
Gamma_average_order=log(Gamma_average_save(1:end-1)./Gamma_average_save(2:end))./log(stepsize(1:end-1)./stepsize(2:end))