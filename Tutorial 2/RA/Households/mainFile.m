%% MIT MATRICES IN GROWTH MODELS
% Christian Wolf
% this version: 03/21/2018

%% HOUSEKEEPING

clc
clear all
close all

%% SHOCK SEQUENCE

global z T

T = 250;

%% PARAMETERS AND STEADY STATE

global beta gamma varphi chi alpha delta rho_z sigma_z ...
    C_SS L_SS W_SS R_SS Y_SS K_SS I_SS A_SS Z_SS

% parameters

beta = 0.99;
gamma = 1;
varphi = 0.5;
alpha = 1/3;
delta = 0.025;
rho_z = 0.9;
sigma_z = 0.1;

% steady state

Z_SS = 1;
L_SS = 1/3;
R_SS = 1/(beta);
K_SS = (alpha * Z_SS/(R_SS - (1-delta)))^(1/(1-alpha)) * L_SS;
W_SS = (1-alpha) * Z_SS * K_SS^alpha * L_SS^(-alpha);
Y_SS = Z_SS * K_SS^alpha * L_SS^(1-alpha);
I_SS = (K_SS - (1-delta) * K_SS);
C_SS = Y_SS - I_SS;
A_SS = K_SS;

chi = W_SS * C_SS^(-gamma) * L_SS^(-1/varphi);

% shock sequences

e_z = zeros(T,1); e_z(1) = 1;

z = zeros(T,1);
for t = 1:T
    if t == 1
        z(t) = sigma_z * e_z(1);
    else
        z(t) = rho_z * z(t-1);
    end
end

%% COMPUTE EQUILIBRIUM

k_guess  = zeros(T,1);
l_guess  = zeros(T,1);
guess    = [k_guess;l_guess];
guess    = myAD(guess);
mktcl_0  = mktcl_fun(zeros(2*T,1));
mktcl    = mktcl_fun(guess);

A = getderivs(mktcl);
A = full(A);
sol = -A^(-1) * mktcl_0;

%% CHECK THAT EVERYTHING IS FINE

%----------------------------------------------------------------
% Collect Results
%----------------------------------------------------------------

k_seq = sol(1:T,1);
l_seq = sol(T+1:2*T,1);

%----------------------------------------------------------------
% Get Prices
%----------------------------------------------------------------

% wages

w_seq = z + alpha * [0;k_seq(1:end-1)] - alpha * l_seq;

% real rate

r_seq = (alpha * Z_SS * K_SS^(alpha-1) * L_SS^(1-alpha))/R_SS * (z + (alpha-1) * [0;k_seq(1:end-1)] + (1-alpha) * l_seq);

%----------------------------------------------------------------
% Solve for Household Behavior
%----------------------------------------------------------------

A = zeros(2*T,2*T);
for t = 1:T
    A(t,t) = 1;
    if t < T
    A(t,t+1) = -1;
    end
    
    A(T+t,t) = C_SS;
    if t > 1
    A(T+t,T+t-1) = -R_SS * A_SS;
    end
    A(T+t,T+t) = A_SS;
end

b_vec = [-1/gamma * [r_seq(2:end);0]; W_SS * L_SS * (w_seq + l_seq) + R_SS * A_SS * r_seq];

sol_seq = A^(-1) * b_vec;
c_seq   = sol_seq(1:T);
a_seq   = sol_seq(T+1:2*T);
k_seq   = a_seq;

y_seq   = z + alpha * [0;k_seq(1:end-1)] + (1-alpha) * l_seq;