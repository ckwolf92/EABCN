function supply = supply_fun(c);

global beta gamma delta alpha nu kappa R_SS L_SS Z_SS K_SS W_SS Y_SS C_SS I_SS Lambda_SS chi T z

% get prices

lambda = - gamma * c;
w = - lambda;
r = 1/(beta * R_SS) * (lambda - [lambda(2:end);0]);
q = zeros(T,1);

% get returns

pi_I = 1/(1 - (1-alpha) * nu) * z - (1-alpha)*nu/(1 - (1-alpha)*nu) * w;

% get behavior of firms

A = zeros(T,T);
for i = 1:T
    A(i,i) = (kappa * I_SS^2 * 1/delta - beta * (1/beta - (1-delta)) * (alpha * nu - 1 + alpha * nu * (1-alpha) * nu/(1 - (1-alpha) * nu)) + beta * (1-delta) * kappa * I_SS^2 * 1/delta * (1-delta));
end
for i = 2:T
    A(i,i-1) = (-kappa * I_SS^2 * 1/delta * (1-delta));
    A(i-1,i) = (- beta * (1-delta) * kappa * I_SS^2 * 1/delta);
end

b = lambda * (-1) + [lambda(2:end);0] * (beta * (1/beta - (1-delta)) + beta * (1-delta)) ...
    + [z(2:end);0] * (beta * (1/beta - (1-delta)) * (1 + (1-alpha)*nu/(1-(1-alpha)*nu))) ...
    + [w(2:end);0] * (- beta * (1/beta - (1-delta)) * (1-alpha)*nu/(1-(1-alpha)*nu));

k = A^(-1) * b;
k = [0;k(1:end-1)];
l = 1/(1 - (1-alpha) * nu) * (z - w) + alpha * nu/(1-(1-alpha)*nu) * k;
y = z + alpha * nu * k + (1-alpha) * nu * l;
i = 1/delta * [k(2:end) - (1-delta) * k(1:end-1);0];

% get net output supply

supply = Y_SS/C_SS * y - I_SS/C_SS * i;