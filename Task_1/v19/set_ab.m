clc; clear;
alpha=(1:1:10)';
x=0.5:0.01:1.5;
y= sin(log(sqrt( 1 + ( 1 - exp( -abs(x) ) ).^2))) .^2 + alpha*log(x);

plot(x, y, LineWidth=1.5);
legend("1", "2", "3", "4", "5", "6", "7", "8", "9", "10");
grid on;

