clear all; clc; close all;
Data = csvread("calc3_Kaplina_u_vs_x.csv", 0, 0);
x = Data(:, 1);
u = Data(:, 2);

plot(x, u);