clear all; clc; close all;
Data = csvread("u.csv", 0, 0);
x = Data(:, 1);
u = Data(:, 2);

plot(x, u);