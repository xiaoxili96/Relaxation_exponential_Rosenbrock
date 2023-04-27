stepsize=[1/100 1/200 1/400 1/800 1/1600]; 

TIME=[];  ERR=[];  ENERGY_ERR=[];
for k=1:size(stepsize,2)
    tau=stepsize(k); 
    [time,err,energy_err]=AVF(tau);
    TIME=[TIME time];
    ERR=[ERR err]; 
    ENERGY_ERR=[ENERGY_ERR energy_err];
end

ERR_AVF=ERR;
Err_Order_AVF=log(ERR(1:end-1)./ERR(2:end))./log(stepsize(1:end-1)./stepsize(2:end));
ENERGY_ERR_AVF=ENERGY_ERR;
TIME_AVF=TIME;
save('data_AVF.mat','ERR_AVF','Err_Order_AVF','ENERGY_ERR_AVF','TIME_AVF');