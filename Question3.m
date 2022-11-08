function [y] = Question3(xo)
x=[0,1.2,2.4,3.6,4.8,6];   % I define our x values.
y=[10.0000,6.7113,0.0251,-3.9177,-2.4043,1.5065];  % I define our y values.
%If we want take input from keybord, we can use:
%"  xo=input('xo: '); " code line.
n=length(x);    % Number of terms of x (=y).
a(1)=y(1);      % Set index 1 to y's index 1.
%Newton divided difference calculate algorithm:
for i=1:n-1      
    L(i,1)=(y(i+1)-y(i))/(x(i+1)-x(i));  % We will use for the following formula.
end
for j=2:n-1    % We start with 2 because of undefinedness.
    for i=1:n-j
        L(i,j)=(L(i+1,j-1)-L(i,j-1))/(x(i+j)-x(i));
    end
end
for j=2:n      % We start with 2 because of undefinedness.
    a(j)=L(1,j-1);  % We advance the values by assigning L to the jth value of a.    
end
y=a(1);
xn=1;        
for k=2:n
    xn=xn*(xo-x(k-1)); 
    y=y+a(k)*xn;      % It calculate y value.  
end
end