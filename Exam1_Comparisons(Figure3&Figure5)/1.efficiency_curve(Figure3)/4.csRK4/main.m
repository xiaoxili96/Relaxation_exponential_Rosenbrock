stepsize=[1/100 1/200 1/400 1/800 1/1600]; 

TIME=[];  ERR=[];  ENERGY_ERR=[];
for k=1:size(stepsize,2)
    tau=stepsize(k); 
    [time,err,energy_err]=csRK4(tau);
    TIME=[TIME time];
    ERR=[ERR err]; 
    ENERGY_ERR=[ENERGY_ERR energy_err];
end

ERR_csRK4=ERR;
Err_Order_csRK4=log(ERR(1:end-1)./ERR(2:end))./log(stepsize(1:end-1)./stepsize(2:end));
ENERGY_ERR_csRK4=ENERGY_ERR;
TIME_csRK4=TIME;
save('data_csRK4.mat','ERR_csRK4','Err_Order_csRK4','ENERGY_ERR_csRK4','TIME_csRK4');
