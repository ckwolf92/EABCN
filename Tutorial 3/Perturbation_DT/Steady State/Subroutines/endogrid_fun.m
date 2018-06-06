function [c,ap] = endogrid_fun(cp,states,L_SS,R_SS,W_SS,Trans_SS,gamma,beta,Emat,fspace);

n_epsi = fspace.n(2);
n_a = fspace.n(1);
grid_a = states(1:n_a,1);
a_lb = min(grid_a);

c_endo = (beta * R_SS * (Emat * cp.^(-gamma))).^(-1/gamma);
a_endo = 1/R_SS * (states(:,1) + c_endo - W_SS * L_SS * states(:,2) - Trans_SS * states(:,2));

ap = NaN * a_endo;
for i_e = 1:n_epsi
    V_temp = interpOneD(grid_a,a_endo((i_e-1)*n_a+1:i_e*n_a));
    ap((i_e-1)*n_a+1:i_e*n_a) = V_temp * grid_a; % on a_endo I map into grid_a, so what do I map into from grid_a?
end
ap = max(a_lb,ap);
c = W_SS * L_SS * states(:,2) + R_SS * states(:,1) + Trans_SS * states(:,2) - ap;

end