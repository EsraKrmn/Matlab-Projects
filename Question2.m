function [y] = Question2(xo)
x=[0,1.2,2.4,3.6,4.8,6];   % I define our x values.
fy=[10.0000,6.7113,0.0251,-3.9177,-2.4043,1.5065];  % I define our y values.

% I create a matrix. We use DIM (in size func. : size(X,DIM)) = 2
% for gain returns a row vector containing the number of rows and columns.

n=size(x,2);             % our x's array size  
L=ones(n,size(xo,2));    % our input's size

% Our inputs and outputs must have equal size each other.
if (size(x,2)~=size(fy,2))
   y=NaN;

% If they are equal:
else
   for i=1:n         % i is start with index 1 and stop at n.
      for j=1:n      % j is start with index 1 and stop at n.
         if (i~=j)   % if i and j are not equal, it returns logical 1.
            % It is calculate Li(x) of lagrange formula:
            L(i,:)=L(i,:).*(xo-x(j))/(x(i)-x(j));  
         end
      end
   end                                                                      
   y=0;
   for i=1:n
      y=y+fy(i)*L(i,:);   % Lagrange interpolation formula
   end
end
end