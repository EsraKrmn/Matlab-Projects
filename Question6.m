x_cnt = 0:0.0001:6;   % We define xcnt vector.               
y = Question2(x_cnt);  

%This is for plotting Lagrange Polinom:
plot(x_cnt,y)

J = besselj(0,x_cnt)*10;  % We define our Yexact.
for i = 0
    J(i+1,:) = besselj(0,x_cnt)*10;
end
abs_error = abs(J-y);  %It calculates absolute error.

hold on %This is for plotting abs_error with Question 2
plot(x_cnt,abs_error)
hold off