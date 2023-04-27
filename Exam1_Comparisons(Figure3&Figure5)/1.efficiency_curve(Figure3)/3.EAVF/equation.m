function result=equation(Zn1_t,W,W_hat,PZ,P1,d,tau,alpha,old,old1)

Un1_t=Zn1_t(d+1:end,1);  new=real(W_hat*Un1_t)-old; 
F=(1/4)*new.^3+old.*new.^2+(1/2)*(3*old.^2.*new+new)+old1;

result=Zn1_t-PZ-tau*(-alpha)*P1*[W*F;zeros(d,1)];
