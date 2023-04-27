function Y=compute_nonlinear(Zmid,d,s,m)

Zmid_reshape=reshape(Zmid,d,s);  Xmid=Zmid_reshape(2*m+1:4*m,:);
a=Xmid(2:m,:)-Xmid(1:m-1,:)-Xmid(m+1:2*m-1,:)-Xmid(m:2*m-2,:);  a=a.^3;
y=[-a;zeros(1,s);-a;zeros(1,s)];  y(2:m,:)=y(2:m,:)+a;  y(m:2*m-2,:)=y(m:2*m-2,:)-a;  
xx1=(Xmid(1,:)-Xmid(m+1,:)).^3;  xx2=(Xmid(m,:)+Xmid(2*m,:)).^3;
y(1,:)=y(1,:)+xx1;  y(m+1,:)=y(m+1,:)-xx1;
y(m,:)=y(m,:)+xx2;  y(2*m,:)=y(2*m,:)+xx2;

Y=[-y; zeros(2*m,s)];
Y=Y(:);


