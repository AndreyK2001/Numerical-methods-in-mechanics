clear all; clc; close all;
Data = csvread("n2000.csv", 0, 0);
x = Data(: , 1);
u = Data(:, 2);

plot(x, u);