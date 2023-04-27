load('reference.mat')
stepsize=[1/160 1/320 1/640 1/1280 1/2560]; 

Err=zeros(1,size(stepsize,2));
for k=1:size(stepsize,2)
    tau=stepsize(k); 
    zn=Solver(tau);
    Err(1,k)=max(abs(zn-Zn_c_100000));
end

Err_Order=log(Err(1:end-1)./Err(2:end))./log(stepsize(1:end-1)./stepsize(2:end))
