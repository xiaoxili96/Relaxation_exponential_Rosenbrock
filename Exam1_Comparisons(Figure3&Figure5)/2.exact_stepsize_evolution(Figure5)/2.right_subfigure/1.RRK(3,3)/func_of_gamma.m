function result=func_of_gamma(Un1,Un2,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma)

Un=Un1+gamma*Un2;
result=energy_1 + gamma*energy_2 + gamma^2*energy_3 + h*sum(func_F(Un)) - energy_old;