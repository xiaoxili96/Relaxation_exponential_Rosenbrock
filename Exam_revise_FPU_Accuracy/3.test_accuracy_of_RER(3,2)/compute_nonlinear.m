function Y=compute_nonlinear(X,m)

a=X(2:m)-X(1:m-1)-X(m+1:2*m-1)-X(m:2*m-2);  a=a.^3;
y=[-a;0;-a;0];  y(2:m)=y(2:m)+a;  y(m:2*m-2)=y(m:2*m-2)-a;  
xx1=(X(1)-X(m+1)).^3;  xx2=(X(m)+X(2*m)).^3;
y(1)=y(1)+xx1;  y(m+1)=y(m+1)-xx1;
y(m)=y(m)+xx2;  y(2*m)=y(2*m)+xx2;

Y=[-y; zeros(2*m,1)];
Y=Y(:);


