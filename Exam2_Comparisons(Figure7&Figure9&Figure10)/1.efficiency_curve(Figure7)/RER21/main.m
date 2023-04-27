stepsize=[1/2 1/4 1/8 1/16 1/32];

TIME=[];  ERR=[];  ENERGY_ERR=[];
for k=1:size(stepsize,2)
    tau=stepsize(k); 
    [time,err,energy_err]=rexprb2s1(tau);
    TIME=[TIME time];
    ERR=[ERR err]; 
    ENERGY_ERR=[ENERGY_ERR energy_err];
end

Err_Order=log(ERR(1:end-1)./ERR(2:end))./log(stepsize(1:end-1)./stepsize(2:end))
