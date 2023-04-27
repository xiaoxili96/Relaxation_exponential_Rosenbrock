function Der=compute_nonlinear_der(X,m)

o2m=zeros(2*m,2*m);
a=X(2:m)-X(1:m-1)-X(m+1:2*m-1)-X(m:2*m-2);  a=3*a.^2;
b=3*(X(1)-X(m+1)).^2;  c=3*(X(m)+X(2*m)).^2;

A_der1=-diag(a)+diag(a(1:m-2),1);
A_der2=-diag(a)-diag(a(1:m-2),1);
A_der=[A_der1 A_der2 zeros(m-1,2)];  A_der(m-1,m)=a(m-1);  A_der(m-1,2*m-1)=-a(m-1);
F_der=[-A_der; zeros(1,2*m);-A_der; zeros(1,2*m)];
F_der(2:m,:)=F_der(2:m,:)+A_der;  F_der(m:2*m-2,:)=F_der(m:2*m-2,:)-A_der;
F_der(1,1)=F_der(1,1)+b;  F_der(1,m+1)=F_der(1,m+1)-b;
F_der(m,m)=F_der(m,m)+c;  F_der(m,2*m)=F_der(m,2*m)+c; 
F_der(m+1,1)=F_der(m+1,1)-b;  F_der(m+1,m+1)=F_der(m+1,m+1)+b;
F_der(2*m,m)=F_der(2*m,m)+c;  F_der(2*m,2*m)=F_der(2*m,2*m)+c; 
F_der1=[o2m F_der; o2m o2m];

Der=-F_der1;  


