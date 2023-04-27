function result=func_of_energy(Vn_t,Un_t,W_hat,K,area,h,func_F)

Un=real(W_hat*Un_t);
result=0.5*area*real(Vn_t'*Vn_t) ... 
      -0.5*area*real(Un_t'*K*Un_t) ... 
      +h*sum(func_F(Un));