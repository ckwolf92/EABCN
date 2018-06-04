function [c_opt,ap_opt,mutilde_upd] = EGP_fun(mutilde_opt,r_hat_SS,W_SS,L_SS,Trans_SS,beta_hat,gamma,tau_L,states,grid_a,n_yT,Emat_yP);

% preparations

n_a = size(grid_a,2);
n_y = size(states,1)/n_a;

% re-scale

mutilde_opt_rep = repmat(mutilde_opt,n_yT,1);

% consumption and savings on endogenous grid
    
c_endo = (beta_hat * (1 + r_hat_SS) * mutilde_opt_rep).^(-1/gamma);
a_endo = 1/(1+r_hat_SS) * (states(:,1) + c_endo - (1-tau_L) * W_SS * L_SS * states(:,2) - Trans_SS);

% interpolate

a_endo = reshape(a_endo,n_a,n_y);
V_temp = interpOneD_vec(repmat(grid_a',1,n_y),a_endo);
ap_temp = V_temp * repmat(grid_a',n_y,1);

% deal with constraint

constr = (ap_temp < 0);
ap_opt = constr .* 0 + (1-constr) .* ap_temp;

% get final consumption

c_opt = (1 + r_hat_SS) * states(:,1) + (1-tau_L) * W_SS * L_SS * states(:,2) + Trans_SS - ap_opt;
u_c_opt = c_opt.^(-gamma);

% get mutilde on exogenous grid

ap_opt_temp = reshape(ap_opt,n_a,n_y);
V_temp = interpOneD_vec(ap_opt_temp,repmat(grid_a',1,n_y));
mutilde_on_ap = V_temp * mutilde_opt_rep;

% approximate multiplier

mu_opt = u_c_opt - beta_hat * (1 + r_hat_SS) * mutilde_on_ap;
mu_opt = constr .* mu_opt;
mu_opt = max(0,mu_opt);

% get updated mu_tilde

mutilde_upd = Emat_yP * (c_opt.^(-gamma));

end