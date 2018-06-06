function obj_VFI = obj_VFI_fun_t(c,states_VFI,L,R,W,Trans,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);

ap    = W .* L .* states_VFI(:,2) + R * states_VFI(:,1) + Trans * states_VFI(:,2) - c;
Phi_Ap_VFI = splibas(grid_a_VFI_0,0,spliorder(1),ap);

if gamma == 1
    payoff_today = log(c);
else
    payoff_today = (c.^(1-gamma) - 1)./(1-gamma);
end
payoff_tomorrow = dprod(Phi_epsi_VFI,Phi_Ap_VFI) * h2;

obj_VFI = payoff_today + beta * payoff_tomorrow;

end