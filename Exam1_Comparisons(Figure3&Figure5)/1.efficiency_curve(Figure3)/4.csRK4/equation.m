function result=equation(ZZ1,W,W_hat,A,d,tau,alpha,coe1,coe2,coe3,coe4,Y0,YY0,YY02,YY03)

Y1=ZZ1(1:2*d,1);  Y2=ZZ1(2*d+1:end,1);
YY1=real(W_hat*Y1(d+1:end,1));  YY2=real(W_hat*Y2(d+1:end,1));  
YY12=YY1.^2;  YY13=YY1.^3;  YY22=YY2.^2;  YY23=YY2.^3;
X1=YY03;  
X2=YY02.*YY1;  
X3=YY02.*YY2;  
X4=YY0.*YY12;  
X5=YY0.*YY1.*YY2;
X6=YY0.*YY22;  
X7=YY0;  
X8=YY13;  
X9=YY12.*YY2;  
X10=YY1.*YY22;  
X11=YY1;  
X12=YY23;  
X13=YY2; 

F1=coe3(1)*X1+coe3(2)*X2+coe3(3)*X3+coe3(4)*X4+coe3(5)*X5+coe3(6)*X6 ...
  +coe3(7)*X7+coe3(8)*X8+coe3(9)*X9+coe3(10)*X10+coe3(11)*X11+coe3(12)*X12 ...
  +coe3(13)*X13;
F2=coe4(1)*X1+coe4(2)*X2+coe4(3)*X3+coe4(4)*X4+coe4(5)*X5+coe4(6)*X6 ...
  +coe4(7)*X7+coe4(8)*X8+coe4(9)*X9+coe4(10)*X10+coe4(11)*X11+coe4(12)*X12 ...
  +coe4(13)*X13;

result=ZZ1-[Y0;Y0] ... 
      -tau*[A*(coe1(1)*Y0+coe1(2)*Y1+coe1(3)*Y2);A*(coe2(1)*Y0+coe2(2)*Y1+coe2(3)*Y2)] ...
      -tau*(-alpha)*[W*F1;zeros(d,1);W*F2;zeros(d,1)];