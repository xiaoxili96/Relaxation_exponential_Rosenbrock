function result=func_of_energy(energy_1,Un,h,func_F)

result=energy_1 + h*h*sum(func_F(Un));