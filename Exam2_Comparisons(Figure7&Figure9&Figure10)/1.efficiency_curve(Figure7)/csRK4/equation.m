function result=equation(ZZ1,Zn,Un,Yn,func_f,A,d,tau,coe1,coe2,coe3,coe4,coe5,coe6,GW)

Zm=ZZ1(1:2*d,1);  Zn1=ZZ1(2*d+1:end,1);
Um=Zm(d+1:end,1);  Un1=Zn1(d+1:end,1);

Ym=A*Zm;  Yn1=A*Zn1;
L1=coe5(1)*Yn+coe5(2)*Ym+coe5(3)*Yn1;
L2=coe6(1)*Yn+coe6(2)*Ym+coe6(3)*Yn1;

FF1=func_f(coe1(1)*Un+coe2(1)*Um+coe3(1)*Un1);
FF2=func_f(coe1(2)*Un+coe2(2)*Um+coe3(2)*Un1);
FF3=func_f(coe1(3)*Un+coe2(3)*Um+coe3(3)*Un1);
FF4=func_f(coe1(4)*Un+coe2(4)*Um+coe3(4)*Un1);
FF5=func_f(coe1(5)*Un+coe2(5)*Um+coe3(5)*Un1);
F1=coe4(1)*FF1+coe4(2)*FF2+coe4(3)*FF3+coe4(4)*FF4+coe4(5)*FF5;  F1=(-1)*F1;
F2=GW(1)*FF1+GW(2)*FF2+GW(3)*FF3+GW(4)*FF4+GW(5)*FF5;  F2=(-1)*F2;


result=ZZ1-[Zn;Zn]-tau*[L1;L2]-tau*[F1;zeros(d,1);F2;zeros(d,1)];