function Y_der=compute_nonlinear_der(Zmid,d,s,m)

o2m=zeros(2*m,2*m);
Zmid_reshape=reshape(Zmid,d,s);  Xmid=Zmid_reshape(2*m+1:4*m,:);
a=Xmid(2:m,:)-Xmid(1:m-1,:)-Xmid(m+1:2*m-1,:)-Xmid(m:2*m-2,:);  a=3*a.^2;
b=3*(Xmid(1,:)-Xmid(m+1,:)).^2;  c=3*(Xmid(m,:)+Xmid(2*m,:)).^2;

g=a(:,1);  e=b(1);  f=c(1);
A_der1=-diag(g)+diag(g(1:m-2),1);
A_der2=-diag(g)-diag(g(1:m-2),1);
A_der=[A_der1 A_der2 zeros(m-1,2)];  A_der(m-1,m)=g(m-1);  A_der(m-1,2*m-1)=-g(m-1);
F_der=[-A_der; zeros(1,2*m);-A_der; zeros(1,2*m)];
F_der(2:m,:)=F_der(2:m,:)+A_der;  F_der(m:2*m-2,:)=F_der(m:2*m-2,:)-A_der;
F_der(1,1)=F_der(1,1)+e;  F_der(1,m+1)=F_der(1,m+1)-e;
F_der(m,m)=F_der(m,m)+f;  F_der(m,2*m)=F_der(m,2*m)+f; 
F_der(m+1,1)=F_der(m+1,1)-e;  F_der(m+1,m+1)=F_der(m+1,m+1)+e;
F_der(2*m,m)=F_der(2*m,m)+f;  F_der(2*m,2*m)=F_der(2*m,2*m)+f; 
F_der1=[o2m F_der; o2m o2m];

g=a(:,2);  e=b(2);  f=c(2);
A_der1=-diag(g)+diag(g(1:m-2),1);
A_der2=-diag(g)-diag(g(1:m-2),1);
A_der=[A_der1 A_der2 zeros(m-1,2)];  A_der(m-1,m)=g(m-1);  A_der(m-1,2*m-1)=-g(m-1);
F_der=[-A_der; zeros(1,2*m);-A_der; zeros(1,2*m)];
F_der(2:m,:)=F_der(2:m,:)+A_der;  F_der(m:2*m-2,:)=F_der(m:2*m-2,:)-A_der;
F_der(1,1)=F_der(1,1)+e;  F_der(1,m+1)=F_der(1,m+1)-e;
F_der(m,m)=F_der(m,m)+f;  F_der(m,2*m)=F_der(m,2*m)+f; 
F_der(m+1,1)=F_der(m+1,1)-e;  F_der(m+1,m+1)=F_der(m+1,m+1)+e;
F_der(2*m,m)=F_der(2*m,m)+f;  F_der(2*m,2*m)=F_der(2*m,2*m)+f; 
F_der2=[o2m F_der; o2m o2m];

% g=a(:,3);  e=b(3);  f=c(3);
% A_der1=-diag(g)+diag(g(1:m-2),1);
% A_der2=-diag(g)-diag(g(1:m-2),1);
% A_der=[A_der1 A_der2 zeros(m-1,2)];  A_der(m-1,m)=g(m-1);  A_der(m-1,2*m-1)=-g(m-1);
% F_der=[-A_der; zeros(1,2*m);-A_der; zeros(1,2*m)];
% F_der(2:m,:)=F_der(2:m,:)+A_der;  F_der(m:2*m-2,:)=F_der(m:2*m-2,:)-A_der;
% F_der(1,1)=F_der(1,1)+e;  F_der(1,m+1)=F_der(1,m+1)-e;
% F_der(m,m)=F_der(m,m)+f;  F_der(m,2*m)=F_der(m,2*m)+f; 
% F_der(m+1,1)=F_der(m+1,1)-e;  F_der(m+1,m+1)=F_der(m+1,m+1)+e;
% F_der(2*m,m)=F_der(2*m,m)+f;  F_der(2*m,2*m)=F_der(2*m,2*m)+f; 
% F_der3=[o2m F_der; o2m o2m];
 
% g=a(:,4);  e=b(4);  f=c(4);
% A_der1=-diag(g)+diag(g(1:m-2),1);
% A_der2=-diag(g)-diag(g(1:m-2),1);
% A_der=[A_der1 A_der2 zeros(m-1,2)];  A_der(m-1,m)=g(m-1);  A_der(m-1,2*m-1)=-g(m-1);
% F_der=[-A_der; zeros(1,2*m);-A_der; zeros(1,2*m)];
% F_der(2:m,:)=F_der(2:m,:)+A_der;  F_der(m:2*m-2,:)=F_der(m:2*m-2,:)-A_der;
% F_der(1,1)=F_der(1,1)+e;  F_der(1,m+1)=F_der(1,m+1)-e;
% F_der(m,m)=F_der(m,m)+f;  F_der(m,2*m)=F_der(m,2*m)+f; 
% F_der(m+1,1)=F_der(m+1,1)-e;  F_der(m+1,m+1)=F_der(m+1,m+1)+e;
% F_der(2*m,m)=F_der(2*m,m)+f;  F_der(2*m,2*m)=F_der(2*m,2*m)+f; 
% F_der4=[o2m F_der; o2m o2m];

Y_der=zeros(d*s,d*s);
Y_der(1:d,1:d)=-F_der1;  
Y_der(d+1:2*d,d+1:2*d)=-F_der2;  
% Y_der(2*d+1:3*d,2*d+1:3*d)=-F_der3;
% Y_der(3*d+1:4*d,3*d+1:4*d)=-F_der4;


