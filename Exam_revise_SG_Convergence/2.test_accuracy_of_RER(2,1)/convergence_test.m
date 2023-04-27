stepsize=[1/512 1/1024 1/2048];

V_save=[];  U_save=[];  Gamma_average_save=[];
for k=1:size(stepsize,2)
    tau=stepsize(1,k);
    [Vn,Un,gamma_average]=rexprb2s1(tau);
    V_save=[V_save Vn];  U_save=[U_save Un];  Gamma_average_save=[Gamma_average_save gamma_average];
end

err=[];
for k=1:(size(U_save,2)-1)
    err=[err max([abs(V_save(:,k)-V_save(:,end));abs(U_save(:,k)-U_save(:,end))])];
end

err 
Err_order=log(err(1:end-1)./err(2:end))./log(stepsize(1:end-2)./stepsize(2:end-1))
Gamma_average_save
Gamma_average_order=log(Gamma_average_save(1:end-1)./Gamma_average_save(2:end))./log(stepsize(1:end-1)./stepsize(2:end))