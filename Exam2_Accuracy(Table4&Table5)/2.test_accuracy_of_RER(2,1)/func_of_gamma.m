function result=func_of_gamma(Vn,Un,K,h,func_F,Vn_update,Un_update,energy_old,gamma)

result=func_of_energy(Vn+gamma*Vn_update,Un+gamma*Un_update,K,h,func_F)-energy_old;