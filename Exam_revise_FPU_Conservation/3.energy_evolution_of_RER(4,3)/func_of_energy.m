function result=func_of_energy(yn,xn,AA,m)

an=xn(2:m)-xn(1:m-1)-xn(m+1:2*m-1)-xn(m:2*m-2);  
result=0.5*(yn'*yn)+0.5*(xn'*AA*xn)+0.25*(sum(an.^4)+(xn(1)-xn(m+1))^4+(xn(m)+xn(2*m))^4);