function result=func_of_gamma(Vn_t,Un_t,W_hat,K,area,h,func_F,Vn_t_update,Un_t_update,energy_old,gamma)

result=func_of_energy(Vn_t+gamma*Vn_t_update,Un_t+gamma*Un_t_update,W_hat,K,area,h,func_F)-energy_old;