function vResidual = equilibriumConditions_reduced_AD(vars)

%% PREPARATIONS

%----------------------------------------------------------------
% Global Variables
%----------------------------------------------------------------

% aggregate parameters

global beta gamma varphi chi epsilon_p theta_p rho_tr phi_pi phi_y ...
    rho_m sigma_m

% household parameters

global a_lb n_epsi grid_epsi epsi_dist Pi_epsi Gamma

% steady-state quantities

global C_SS Y_SS L_SS W_SS A_SS Trans_SS MC_SS Pi_SS R_n_SS ...
    Wedge_EE_SS c_opt_SS ap_opt_SS mu_tilde_SS lambda_SS lambda_orig_SS

% other quantities

global aux grid_a grid_a_0 spliorder tol_GS_choice tol_GS_VFI tol_GS_dist ...
    states Phi_epsi Phi Emat fspaceerga fspace ...
    n_a_VFI n_a_dist n_s_VFI n_s_dist a_min a_max n_a_VFI_reduced n_s_VFI_reduced n_s_dist_reduced

% system reduction

global from_spline to_spline reduce_VF reduce_dist lambda_epsi_SS lambda_a_SS

% system size

global n_shocks nVars nEErrors vars_SS

%----------------------------------------------------------------
% Unpack Variables 
%----------------------------------------------------------------

% current values

mu_tilde = vars_SS(1:n_s_VFI) + from_spline * vars(1:n_s_VFI_reduced);

if reduce_dist == 0

    lambda = vars_SS(n_s_VFI+1:n_s_VFI+n_s_dist) + vars(n_s_VFI_reduced+1:n_s_VFI_reduced+n_s_dist_reduced);
    
elseif reduce_dist == 1
    
    lambda_epsi = lambda_epsi_SS;
    lambda_a    = vars_SS(n_s_VFI+1:n_s_VFI+n_a_dist) + vars(n_s_VFI_reduced+1:n_s_VFI_reduced+n_a_dist);
    
    lambda = copula_reduction(lambda_epsi,lambda_a,lambda_SS);
end

c_agg     = vars(n_s_VFI_reduced+n_s_dist_reduced+1);
y_agg     = vars(n_s_VFI_reduced+n_s_dist_reduced+2);
l_agg     = vars(n_s_VFI_reduced+n_s_dist_reduced+3);
w_agg     = vars(n_s_VFI_reduced+n_s_dist_reduced+4);
a_agg     = vars(n_s_VFI_reduced+n_s_dist_reduced+5);
mc_agg    = vars(n_s_VFI_reduced+n_s_dist_reduced+6);
pi_agg    = vars(n_s_VFI_reduced+n_s_dist_reduced+7);
r_n_agg   = vars(n_s_VFI_reduced+n_s_dist_reduced+8);
trans_agg = vars(n_s_VFI_reduced+n_s_dist_reduced+9);
u_m_agg   = vars(n_s_VFI_reduced+n_s_dist_reduced+10);

% future values

mu_tilde_p1 = vars_SS(1:n_s_VFI) + from_spline * vars(nVars+1:nVars+n_s_VFI_reduced);

if reduce_dist == 0

    lambda_p1 = vars_SS(n_s_VFI+1:n_s_VFI+n_s_dist) + vars(nVars+n_s_VFI_reduced+1:nVars+n_s_VFI_reduced+n_s_dist_reduced);
    
elseif reduce_dist == 1
    
    lambda_epsi_p1 = lambda_epsi_SS;
    lambda_a_p1    = vars_SS(n_s_VFI+1:n_s_VFI+n_a_dist) + vars(nVars+n_s_VFI_reduced+1:nVars+n_s_VFI_reduced+n_a_dist);
    
    lambda_p1 = copula_reduction(lambda_epsi_p1,lambda_a_p1,lambda_SS);
    
end

c_agg_p1     = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+1);
y_agg_p1     = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+2);
l_agg_p1     = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+3);
w_agg_p1     = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+4);
a_agg_p1     = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+5);
mc_agg_p1    = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+6);
pi_agg_p1    = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+7);
r_n_agg_p1   = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+8);
trans_agg_p1 = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+9);
u_m_agg_p1   = vars(nVars+n_s_VFI_reduced+n_s_dist_reduced+10);

% expectational errors

mu_tilde_EError = vars(2*nVars+1:2*nVars+n_s_VFI_reduced,1);
NKPC_EError     = vars(2*nVars+n_s_VFI_reduced+1,1);

m_Shock       = vars(2*nVars+nEErrors+1,1);

% collect aggregates that you will need in levels

W_p1     = W_SS * exp(w_agg_p1);
R_p1     = R_n_SS * exp(r_n_agg)/(Pi_SS * exp(pi_agg_p1));
L_p1     = L_SS * exp(l_agg_p1);
Trans_p1 = Trans_SS * exp(trans_agg_p1);

%% COMPUTE RESIDUALS

%----------------------------------------------------------------
% Household Block
%----------------------------------------------------------------

[c_opt_t,ap_opt_t] = endogrid_fun_t(mu_tilde_p1,states,L_p1,R_p1,W_p1,Trans_p1,gamma,beta,fspace);
mu_tilde_out = Emat * c_opt_t.^(-gamma) + from_spline * mu_tilde_EError;

% update distribution

QZ = sparse(kron(Pi_epsi,ones(n_a_dist,1)));

ap_opt_t   = max(min(ap_opt_t,a_max),a_min);
QAt        = funbas_AD(fspaceerga,ap_opt_t);
        
Qt  = dprod_AD(QZ,QAt);

if reduce_dist == 0
    lambda_p1_out     = Qt' * lambda;    
elseif reduce_dist == 1    
    lambda_p1_out     = Qt' * reshape(lambda',n_a_dist*n_epsi,1);    
end
lambda_p1_out_aux = permute(reshape(lambda_p1_out,[n_a_dist,n_epsi]),[2 1]);

% compute aggregates

c_opt_t  = permute(reshape(c_opt_t,[n_a_dist,n_epsi]),[2 1]);
ap_opt_t = permute(reshape(ap_opt_t,[n_a_dist,n_epsi]),[2 1]);

A_grid = ap_opt_t;
C_grid = c_opt_t;

if reduce_dist == 0
    lambda_aux = permute(reshape(lambda,[n_a_dist,n_epsi]),[2 1]);
elseif reduce_dist == 1    
    lambda_aux = lambda;   
end

A_p1_out  = sum(sum(A_grid(:,:).*lambda_aux(:,:)));
C_p1_out  = sum(sum(C_grid(:,:).*lambda_aux(:,:)));

% collect residuals

mu_tilde_Residual = to_spline * (mu_tilde_out - mu_tilde);

if reduce_dist == 0
    
    lambda_Residual = lambda_p1_out - lambda_p1;
    
elseif reduce_dist == 1

    lambda_a_p1_out    = sum(lambda_p1_out_aux,1)';
    lambda_Residual    = lambda_a_p1_out - lambda_a_p1;

end

BC_Residual = A_p1_out - a_agg_p1;
C_Residual  = log(C_p1_out./C_SS) - c_agg_p1;

%----------------------------------------------------------------
% Rest of Economy
%----------------------------------------------------------------

% labor equation

L_Residual = chi * (L_SS * exp(l_agg_p1))^(1/varphi) - W_SS * exp(w_agg_p1) * (C_SS * exp(c_agg_p1))^(-gamma);

% household aggregation

Wealth_Residual = A_SS + a_agg_p1;

% intermediate goods producer

Output_Residual = Y_SS * exp(y_agg_p1) - L_SS * exp(l_agg_p1);

% retailers

MC_Residual = MC_SS * exp(mc_agg_p1) - W_SS * exp(w_agg_p1);
NKPC_Residual = Pi_SS * exp(pi_agg) * (Pi_SS * exp(pi_agg) - 1) - ((1-epsilon_p)/theta_p + (epsilon_p)/theta_p * MC_SS * exp(mc_agg) ...
    - 0.5 * (1-epsilon_p) * (Pi_SS * exp(pi_agg) - 1)^2 + beta * (Pi_SS * exp(pi_agg_p1) - 1) * Pi_SS * exp(pi_agg_p1) * exp(y_agg_p1 - y_agg) + NKPC_EError);
Trans_Residual = Trans_SS * exp(trans_agg_p1) - ((1 - MC_SS * exp(mc_agg_p1)) * Y_SS * exp(y_agg_p1));

% government

TR_Residual = r_n_agg_p1 - (rho_tr * r_n_agg + (1-rho_tr) * (phi_pi * pi_agg_p1 + phi_y * y_agg_p1) + u_m_agg_p1);

% aggregation

mktcl_Residual = Y_SS * exp(y_agg_p1) - C_SS * exp(c_agg_p1);

% shocks

m_Residual = u_m_agg_p1 - rho_m * u_m_agg - sigma_m * m_Shock;

% collect all results

vResidual = [mu_tilde_Residual;lambda_Residual;...
    C_Residual;L_Residual;...
    BC_Residual;Wealth_Residual;...
    Output_Residual;...
    MC_Residual;NKPC_Residual;Trans_Residual;...
    TR_Residual;...
    m_Residual];

end