function result=func_of_gamma(Un,Un_Update,energy_1,energy_2,energy_3,h,func_F,energy_old,gamma)

Un1=Un+gamma*Un_Update;
result=energy_1 + gamma*energy_2 + gamma^2*energy_3 + h*h*sum(func_F(Un1)) - energy_old;