function [x,k]=Jacobi(A,b,epsilon)
x1=0;
x2=0;
x3=0;
k=0;
x1e = inf;
x2e = inf;
x3e = inf;

if size(A) ~= 3
    error('The given matrix is not 3 by 3')
    
elseif (abs(A(1,2))+abs(A(1,3))>abs(A(1,1)) || abs((A(2,1))+abs(A(2,3)))>abs(A(2,2)) || abs((A(3,1))+abs(A(3,2)))>abs(A(3,3)))
    error('The given matrix is not diagonally dominant')

else    
while (x1e > epsilon) && (x2e > epsilon) && (x3e > epsilon)
    
x1k=1/A(1,1)*(b(1)-A(1,2)*x2-A(1,3)*x3);
x2k=1/A(2,2)*(b(2)-A(2,1)*x1-A(2,3)*x3);
x3k=1/A(3,3)*(b(3)-A(3,1)*x1-A(3,2)*x2);

x1e = abs(x1k-x1);
x2e = abs(x2k-x2);
x3e = abs(x3k-x3);

x1 = x1k;
x2 = x2k;
x3 = x3k;

k=k+1;
x=[x1,x2,x3];
end
end
end
