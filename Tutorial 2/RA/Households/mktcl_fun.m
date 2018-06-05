function mktcl = mktcl_fun(guess);

global beta gamma varphi chi alpha delta rho_z sigma_z ...
    C_SS L_SS W_SS R_SS Y_SS K_SS I_SS A_SS Z_SS ...
    z T

%----------------------------------------------------------------
% Collect Inputs
%----------------------------------------------------------------

k_in = guess(1:T,1);
l_in = guess(T+1:2*T,1);

%----------------------------------------------------------------
% Get Prices
%----------------------------------------------------------------

% wages

w = z + alpha * [0;k_in(1:end-1)] - alpha * l_in;

% real rate

r = (alpha * Z_SS * K_SS^(alpha-1) * L_SS^(1-alpha))/R_SS * (z + (alpha-1) * [0;k_in(1:end-1)] + (1-alpha) * l_in);

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

b_vec = [-1/gamma * [r(2:end);0]; W_SS * L_SS * (w + l_in) + R_SS * A_SS * r];

sol = A^(-1) * b_vec;
c = sol(1:T);
a = sol(T+1:2*T);

%----------------------------------------------------------------
% Get Distance
%----------------------------------------------------------------

mktcl = [a - k_in; varphi * (w - gamma * c) - l_in];