%% MIT SHOCK SOLUTION FOR SIMPLE GROWTH MODEL
% Christian Wolf
% this version: 09/28/2017

%% HOUSEKEEPING

clc
clear all
close all

%% SHOCK SEQUENCE

global z T

load z
T = size(z,1);

%% PARAMETERS AND STEADY STATE

global beta gamma delta alpha nu kappa R_SS L_SS Z_SS K_SS W_SS Y_SS C_SS I_SS Lambda_SS chi

beta = 0.99;
gamma = 1;
alpha = 1/3;
nu = 0.9;
delta = 0.025;
kappa = 0;

Z_SS    = 1;
R_SS    = 1/beta;
L_SS    = 1/3;

K_SS = ((alpha * nu * Z_SS * L_SS^((1-alpha)*nu))/(1/beta - (1-delta)))^(1/(1-alpha*nu));
Y_SS = Z_SS * (K_SS^alpha * L_SS^(1-alpha))^nu;

I_SS = delta * K_SS;
C_SS = Y_SS - I_SS;
Lambda_SS = (C_SS)^(-gamma);
W_SS = (1-alpha) * nu * Z_SS * (K_SS)^(alpha * nu) * (L_SS)^((1-alpha) * nu - 1);

chi = W_SS * Lambda_SS;

%% COMPUTE EQUILIBRIUM

c_guess  = zeros(T,1);
guess    = c_guess;
guess    = myAD(guess);
supply_0 = supply_fun(zeros(T,1));
supply   = supply_fun(guess);

A = getderivs(supply);
A = full(A);
c_sol = (eye(T) - A)^(-1) * supply_0;