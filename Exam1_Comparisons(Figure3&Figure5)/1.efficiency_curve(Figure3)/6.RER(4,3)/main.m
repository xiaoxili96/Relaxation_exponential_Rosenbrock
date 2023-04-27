stepsize=[1/200 1/400 1/800 1/1600 1/3200]; 

TIME=[];  ERR=[];  ENERGY_ERR=[];
for k=1:size(stepsize,2)
    tau=stepsize(k); 
    [time,err,energy_err]=rexprb4s3(tau);
    TIME=[TIME time];
    ERR=[ERR err]; 
    ENERGY_ERR=[ENERGY_ERR energy_err];
end

ERR_rexprb4s3=ERR;
Err_Order_rexprb4s3=log(ERR(1:end-1)./ERR(2:end))./log(stepsize(1:end-1)./stepsize(2:end));
ENERGY_ERR_rexprb4s3=ENERGY_ERR;
TIME_rexprb4s3=TIME;
save('data_rexprb4s3.mat','ERR_rexprb4s3','Err_Order_rexprb4s3','ENERGY_ERR_rexprb4s3','TIME_rexprb4s3');
