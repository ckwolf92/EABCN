function obj_VFI = obj_VFI_fun(c,states_VFI,r_hat_SS,Trans_SS,tau_L,gamma,beta_hat,grid_a_VFI_0,spliorder,Phi_y_VFI,h2);

ap         = (1-tau_L) .* states_VFI(:,2) + (1 + r_hat_SS) * states_VFI(:,1) + Trans_SS - c;
Phi_Ap_VFI = splibas(grid_a_VFI_0,0,spliorder(1),ap);

if gamma == 1
    payoff_today = log(c);
else
    payoff_today = (c.^(1-gamma) - 1)./(1-gamma);
end
payoff_tomorrow = dprod(Phi_y_VFI,Phi_Ap_VFI) * h2;

obj_VFI = payoff_today + beta_hat * payoff_tomorrow;

end