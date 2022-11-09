x_cnt = 0:0.0001:6;   % We define xcnt vector.               
y = Question2(x_cnt);  
plot(x_cnt,y)   % We drew the graph for option b. 

hold on    %to add a second line plot without deleting the existing line plot.
y2 = Question3(x_cnt);
plot(x_cnt,y2)   % For option c we plotted the graph in addition to b.
hold off

