clc; clear;
alpha=(1:1:10)';
x=-0.2:0.01:0.8;
y=abs(x .^2 - alpha) - exp(alpha * abs(x));
plot(x, y, LineWidth=2);
legend("1", "2", "3", "4", "5", "6", "7", "8", "9", "10");
grid on;

