clc; clear;
alpha=(1:1:10)';
x=-0.9:0.01:0.1;
y= x.^9 + x.^7 + x.^3 + x - 0.3*alpha*tanh(x);

plot(x, y, LineWidth=1.5);
legend("1", "2", "3", "4", "5", "6", "7", "8", "9", "10");
grid on;

