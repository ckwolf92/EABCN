function [c,ap] = endogrid_fun_t(mu_tilde,states,L,R,W,Trans,gamma,beta,fspace);

n_epsi = fspace.n(2);
n_a = fspace.n(1);
grid_a = states(1:n_a,1);
a_lb = min(grid_a);

c_endo = (beta * R * mu_tilde).^(-1/gamma);
a_endo = 1/R * (states(:,1) + c_endo - W * L * states(:,2) - Trans * states(:,2));

ap = NaN * a_endo;
for i_e = 1:n_epsi
    V_temp = interpOneD(grid_a,a_endo((i_e-1)*n_a+1:i_e*n_a));
    ap((i_e-1)*n_a+1:i_e*n_a) = V_temp * grid_a; % on a_endo I map into grid_a, so what do I map into from grid_a?
end
ap = max(a_lb + 0 * ap,ap);
c = W * L * states(:,2) + R * states(:,1) + Trans * states(:,2) - ap;

end