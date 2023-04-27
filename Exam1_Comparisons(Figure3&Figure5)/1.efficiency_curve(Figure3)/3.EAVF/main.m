stepsize=[1/100 1/200 1/400 1/800 1/1600]; 

TIME=[];  ERR=[];  ENERGY_ERR=[];
for k=1:size(stepsize,2)
    tau=stepsize(k); 
    [time,err,energy_err]=EAVF(tau);
    TIME=[TIME time];
    ERR=[ERR err]; 
    ENERGY_ERR=[ENERGY_ERR energy_err];
end

ERR_EAVF=ERR;
Err_Order_EAVF=log(ERR(1:end-1)./ERR(2:end))./log(stepsize(1:end-1)./stepsize(2:end));
ENERGY_ERR_EAVF=ENERGY_ERR;
TIME_EAVF=TIME;
save('data_EAVF.mat','ERR_EAVF','Err_Order_EAVF','ENERGY_ERR_EAVF','TIME_EAVF');
