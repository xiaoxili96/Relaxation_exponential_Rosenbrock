function result=func_of_gamma(yn,xn,AA,m,yn_update,xn_update,energy_old,gamma)

result=func_of_energy(yn+gamma*yn_update,xn+gamma*xn_update,AA,m)-energy_old;