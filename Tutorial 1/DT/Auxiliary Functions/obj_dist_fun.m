function obj_dist = obj_dist_fun(c,states_dist,r_hat_SS,Trans_SS,tau_L,gamma,beta_hat,grid_a_VFI_0,spliorder,Phi_y_dist,h2);

ap         = (1-tau_L) .* states_dist(:,2) + (1 + r_hat_SS) * states_dist(:,1) + Trans_SS  - c;
Phi_Ap_VFI = splibas(grid_a_VFI_0,0,spliorder(1),ap);

if gamma == 1
    payoff_today = log(c);
else
    payoff_today = (c.^(1-gamma) - 1)./(1-gamma);
end
payoff_tomorrow = dprod(Phi_y_dist,Phi_Ap_VFI) * h2;

obj_dist = payoff_today + beta_hat * payoff_tomorrow;

end