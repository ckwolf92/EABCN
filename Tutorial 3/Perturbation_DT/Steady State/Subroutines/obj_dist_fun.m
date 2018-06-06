function obj_dist = obj_dist_fun(c,states_dist,L_SS,R_SS,W_SS,Trans_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_dist,h2);

ap         = W_SS .* L_SS .* states_dist(:,2) + R_SS * states_dist(:,1) + Trans_SS * states_dist(:,2)  - c;
Phi_Ap_VFI = splibas(grid_a_VFI_0,0,spliorder(1),ap);

if gamma == 1
    payoff_today = log(c);
else
    payoff_today = (c.^(1-gamma) - 1)./(1-gamma);
end
payoff_tomorrow = dprod(Phi_epsi_dist,Phi_Ap_VFI) * h2;

obj_dist = payoff_today + beta * payoff_tomorrow;

end