function obj_dist = obj_dist_fun_t(c,states_dist,L,R,W,Trans,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_dist,h2);

ap    = W .* L .* states_dist(:,2) + R * states_dist(:,1) + Trans * states_dist(:,2) - c;
Phi_Ap_VFI = splibas(grid_a_VFI_0,0,spliorder(1),ap);

if gamma == 1
    payoff_today = log(c);
else
    payoff_today = (c.^(1-gamma) - 1)./(1-gamma);
end
payoff_tomorrow = dprod(Phi_epsi_dist,Phi_Ap_VFI) * h2;

obj_dist = payoff_today + beta * payoff_tomorrow;

end