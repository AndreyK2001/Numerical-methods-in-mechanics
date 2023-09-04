clc; clear;
alpha=(1:1:10)';
x=-6:0.01:6;
y= x + exp(-1 ./ (1 + x.^2)) - alpha + 5;

plot(x, y, LineWidth=1.5);
legend("1", "2", "3", "4", "5", "6", "7", "8", "9", "10");
grid on;

