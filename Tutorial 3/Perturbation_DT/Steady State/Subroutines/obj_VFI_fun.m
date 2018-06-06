function obj_VFI = obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,Trans_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);

ap         = W_SS .* L_SS .* states_VFI(:,2) + R_SS * states_VFI(:,1) + states_VFI(:,2) .* Trans_SS - c;
Phi_Ap_VFI = splibas(grid_a_VFI_0,0,spliorder(1),ap);

if gamma == 1
    payoff_today = log(c);
else
    payoff_today = (c.^(1-gamma) - 1)./(1-gamma);
end
payoff_tomorrow = dprod(Phi_epsi_VFI,Phi_Ap_VFI) * h2;

obj_VFI = payoff_today + beta * payoff_tomorrow;

end