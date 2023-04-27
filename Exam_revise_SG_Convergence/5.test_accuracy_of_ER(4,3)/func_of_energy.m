function result=func_of_energy(Vn,Un,K,h,func_F)

result=0.5*h*h*(Vn'*Vn)-0.5*h*h*(Un'*K*Un)+h*h*sum(func_F(Un));