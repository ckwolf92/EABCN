%% HETEROGENEOUS-HOUSEHOLD MODELS: STEADY-STATE COMPUTATION
% Christian Wolf
% this version: 04/01/2018

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/ckwolf/Dropbox/Research/Miscellaneous/EABCN_Chile/Slides_CKW/codes/Tutorial 3/Perturbation_DT';  
% addpath('/Users/ckwolf/Dropbox/Coding/compecon2011/CEtools')
% addpath('/Users/ckwolf/Dropbox/Coding/compecon2011/CEdemos')
  
addpath([path '/Auxiliary Functions'])
addpath([path '/Steady State/Subroutines'])
cd([ path '/Steady State']);

%% HANK SET-UP

disp('Computing HANK steady state...')

%----------------------------------------------------------------
% Parameters
%----------------------------------------------------------------

% household side

beta      = 0.97;
gamma     = 1;
varphi    = 0.5;
a_lb      = -1;

n_epsi      = 5; % 5 or 7
rho_epsi    = 0.1; % 0.966
sigma_epsi  = 1; % 0.017
p_star_HH   = 0;
[Pi_epsi, grid_epsi] = rouwen_star(rho_epsi, sigma_epsi, n_epsi, p_star_HH);
grid_epsi   = exp(grid_epsi); 
epsi_dist   = ergodicdist(Pi_epsi)'; 
grid_epsi   = grid_epsi/sum(grid_epsi(:).*epsi_dist(:));

Gamma = 0;

% firm side

epsilon_p = 10;
theta_p   = 100;

% government

rho_tr    = 0.75;
phi_pi    = 1.2;
phi_y     = 0;

% shock persistence

rho_m     = 0.3;

% shock volatility

sigma_m     = 0.1;

%----------------------------------------------------------------
% Known Steady-State Objects
%----------------------------------------------------------------

Z_SS     = 1;
L_SS     = 1/3;
Y_SS     = Z_SS * L_SS;
Pi_SS    = 1;
MC_SS    = (epsilon_p-1)/epsilon_p;
W_SS     = MC_SS * Z_SS;
B_SS     = 0;
Trans_SS = (1-W_SS) * Y_SS;

%% NUMERICAL SOLUTION PARAMETERS

%----------------------------------------------------------------
% Asset and Productivity Grids
%----------------------------------------------------------------

% asset grid

n_a_0      = 40; % number of grid points 
a_min      = a_lb;
a_max      = 20;
grid_a_lin = 0; % linear grid for VFI?

% productivity grid

grid_epsi    = grid_epsi;

% spline properties

spliorder = [3 1];

% full grid

n_s_0  = n_a_0 * n_epsi;

%----------------------------------------------------------------
% Tolerance 
%----------------------------------------------------------------

tol_GS_choice = 10^(-6); % tolerance for choice set golden search
tol_VFI       = 10^(-8); % tolerance for unconstrained value function iteration
tol_GS_VFI    = 10^(-6); % tolerance in VFI golden search
tol_GS_dist   = 10^(-6); % tolerance in distribution golden search

%----------------------------------------------------------------
% VFI Solution Set-Up 
%----------------------------------------------------------------

totit_Howard  = 50; % number of Howard improvement steps in VFI
disp_VFIdist  = 0; % display distance in VFI iteration?
do_newton     = 0; % Newton updating?
use_prev      = 0; % use previous VFI approximation? (speed gains, but prevents replicability)

%----------------------------------------------------------------
% Loop Preparations
%----------------------------------------------------------------

R_init      = 1 + (1/beta - 1)/2;
R_guess     = R_init; % convenient for later updating
R_upd       = 0.4; % tempered wage updating
R_tol       = 10^(-5); % tolerance for interest loop
R_it_max    = 50; % maximal number of wage iterations
R_ub        = 1/beta; % upper bound on interest rate (RANK rate)
R_lb        = 1; % lower bound on interest rate

%% MAIN LOOP FOR HANK

for R_it = 1:R_it_max
    
%----------------------------------------------------------------
% Set Distance
%----------------------------------------------------------------

dist = 1;

%----------------------------------------------------------------
% Complete Household Inputs
%----------------------------------------------------------------

R_SS     = R_guess;

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------

% asset grid
                    
if grid_a_lin == 1
    if a_min < 0
        n_a_neg      = round(1/4 * n_a_0);
        n_a_pos      = n_a_0 - n_a_neg;
        grid_a_neg_0 = linspace(a_min,0,n_a_neg);
        grid_a_pos_0 = linspace(0,a_max,n_a_pos);
        grid_a_0     = [grid_a_neg_0,grid_a_pos_0(2:end)];
    else
        grid_a_0 = linspace(a_min,a_max,n_a_0);
    end
else
    coeff_power = 0.9;
    power       = 8;
    if a_min < 0
        n_a_neg      = round(1/4 * n_a_0);
        n_a_pos      = n_a_0 - n_a_neg;
        grid_a_neg_0 = linspace(a_min,0,n_a_neg);
        grid_a_pos_0 = linspace(0,1,n_a_pos);
        grid_a_pos_0 = 0 + (a_max-0)*((1 - coeff_power) * grid_a_pos_0 + coeff_power * (grid_a_pos_0.^power));
        grid_a_0     = [grid_a_neg_0,grid_a_pos_0(2:end)];
    else
        grid_a_0   	 = linspace(0,1,n_a_0);
        grid_a_0	 = a_min + (a_max-a_min)*((1 - coeff_power) * grid_a_0 + coeff_power * (grid_a_0.^power));
    end
end

%----------------------------------------------------------------
% Full Grids
%----------------------------------------------------------------

fspace          = fundef({'spli',grid_a_0,0,spliorder(1)},...
                         {'spli',grid_epsi,0,spliorder(2)});
states_grid  = funnode(fspace);
states       = gridmake(states_grid);

grid_a     = states(states(:,2)==states(1,2),1)';
n_a        = size(grid_a,2);

n_s  = n_a * n_epsi;

Phi_epsi   = splibas(grid_epsi,0,spliorder(2),states(:,2));
Phi_A      = splibas(grid_a_0,0,spliorder(1),states(:,1));
Phi        = dprod(Phi_epsi,Phi_A);
Emat       = kron(Pi_epsi,speye(n_a));

%----------------------------------------------------------------
% Endogenous Gridpoint Iteration
%----------------------------------------------------------------

optset('goldenx','tol',tol_GS_VFI)

if R_it == 1
    mu_tilde  = (W_SS .* L_SS .* states(:,2) + (R_SS - 1) * states(:,1) + states(:,2) .* Trans_SS).^(-gamma);
    mu_tilde  = Emat * mu_tilde;
else
    if use_prev == 0
       mu_tilde  = (W_SS .* L_SS .* states(:,2) + (R_SS - 1) * states(:,1) + states(:,2) .* Trans_SS).^(-gamma);
       mu_tilde  = Emat * mu_tilde;
    end
end

while dist > tol_VFI
    [c_opt,ap_opt] = endogrid_fun_v2(mu_tilde,states,L_SS,R_SS,W_SS,Trans_SS,gamma,beta,fspace);
    dist = norm(Emat * c_opt.^(-gamma) - mu_tilde)/norm(mu_tilde);
    if disp_VFIdist == 1
       disp(dist)
    end
    mu_tilde = Emat * c_opt.^(-gamma);
end

%----------------------------------------------------------------
% Distribution: Final Computation
%----------------------------------------------------------------

ap_opt     = max(min(ap_opt,a_max),a_min);
fspaceerga = fundef({'spli',grid_a,0,1});

QZ = kron(Pi_epsi,ones(n_a,1));

QA = funbas(fspaceerga,ap_opt);
Q  = dprod(QZ,QA);

lambda_SS      = full(ergodicdist(Q));
lambda_orig_SS = lambda_SS;
lambda_SS      = permute(reshape(lambda_SS,[n_a,n_epsi]),[2 1]);

%----------------------------------------------------------------
% Compute Aggregates
%----------------------------------------------------------------

% re-order policy function

c_opt_SS  = permute(reshape(c_opt,[n_a,n_epsi]),[2 1]);
ap_opt_SS = permute(reshape(ap_opt,[n_a,n_epsi]),[2 1]);

% compute other aggregates

A_grid        = NaN(n_epsi,n_a); % asset holdings
for i_epsi = 1:n_epsi
    for i_a = 1:n_a
            A_grid(i_epsi,i_a) = grid_a(i_a);                   
    end
end

A_SS  = sum(sum(A_grid(:,:).*lambda_SS(:,:)));

%----------------------------------------------------------------
% Check Asset Market Clearing
%----------------------------------------------------------------

A_err = A_SS;

%----------------------------------------------------------------
% Update Grid Bound
%----------------------------------------------------------------

% a_max = max(10,grid_a_dist(min(round(1.35 * max(sum((lambda_SS_H > 10^(-8)),2))),n_a_dist)));

% ----------------------------------------------------------------
% Interest Rate Update
% ----------------------------------------------------------------

if A_err < -R_tol
    disp(['Interest Rate: ' num2str(R_guess) ', Interest Rate too low: ' num2str(A_err) ]);
    R_lb = R_guess;
    R_guess = (1-R_upd) * R_guess + R_upd * R_ub;
    Y_guess = Y_SS;
elseif A_err > R_tol
    disp(['Interest Rate: ' num2str(R_guess) ', Interest Rate too high: ' num2str(A_err) ]);
    R_ub = R_guess;
    R_guess = (1-R_upd) * R_guess + R_upd * R_lb;
    Y_guess = Y_SS;
elseif abs(A_err) <= R_tol
    disp(['Steady State Found, Interest Rate = ' num2str(R_guess)]);
    disp(['Output guess: ' num2str(Y_guess) ', Total output produced: ' num2str(Y_SS) ]);
    break
end

end

%% COLLECT RESULTS

% get all aggregates

A_grid        = NaN(n_epsi,n_a); % asset holdings
C_grid        = NaN(n_epsi,n_a); % consumption
for i_epsi = 1:n_epsi
    for i_a = 1:n_a
            A_grid(i_epsi,i_a) = grid_a(i_a);                   
            C_grid(i_epsi,i_a) = c_opt_SS(i_epsi,i_a);
    end
end

A_SS  = sum(sum(A_grid(:,:).*lambda_SS(:,:)));
C_SS  = sum(sum(C_grid(:,:).*lambda_SS(:,:)));

R_n_SS = R_SS;
mu_tilde_SS = mu_tilde;

% get chi from aggregate labor relation

chi  = W_SS * (C_SS)^(-gamma) * 1/(L_SS^(1/varphi));

% compute MPCs

MPC_dist_SS = NaN(n_epsi,n_a);
for i_epsi = 1:n_epsi
    for i_a = 1:n_a
        if i_a < n_a
            MPC_dist_SS(i_epsi,i_a) = 1/R_SS * (c_opt_SS(i_epsi,i_a+1) - c_opt_SS(i_epsi,i_a))/(grid_a(i_a+1) - grid_a(i_a));
        else
            MPC_dist_SS(i_epsi,i_a) = MPC_dist_SS(i_epsi,i_a-1);
        end
    end
end

MPC_SS = sum(MPC_dist_SS(:) .* lambda_SS(:));

%% PLOTS

epsi_plot = round(n_epsi/2);
end_plot  = min(n_a,round(1.1 * max(sum((lambda_SS > 10^(-8)),2))));

%----------------------------------------------------------------
% Stationary Distribution
%----------------------------------------------------------------

figure(1)
set(gca,'FontSize',16)
plot(grid_a(1:end_plot),squeeze(lambda_SS(epsi_plot,1:end_plot)))
xlabel('Assets, a')
ylabel('Mass')
title([ 'Distribution, Type ' num2str(epsi_plot) ])

%----------------------------------------------------------------
% Policy Function
%----------------------------------------------------------------

figure(2)
set(gca,'FontSize',16)
plot(grid_a(1:end_plot),squeeze(c_opt_SS(epsi_plot,1:end_plot)))
xlabel('Assets, a')
title([ 'Consumption Policy Function, Type ' num2str(epsi_plot) ])

%% ANALYZE RANK AND TANK

%----------------------------------------------------------------
% Wedges
%----------------------------------------------------------------

Wedge_EE_SS = 1/(beta * R_SS);

%----------------------------------------------------------------
% RANK
%----------------------------------------------------------------

% identify wedges

% get TANK fraction

frac_H = (MPC_SS - MPC_dist_SS(n_epsi,n_a))/(1 - MPC_dist_SS(n_epsi,n_a));

%% SAVE RESULTS

%----------------------------------------------------------------
% Aggregate Parameters
%----------------------------------------------------------------

save param_agg beta gamma varphi chi epsilon_p theta_p rho_tr phi_pi phi_y ...
    rho_m sigma_m

%----------------------------------------------------------------
% Household Parameters
%----------------------------------------------------------------

save param_households a_lb n_epsi grid_epsi epsi_dist Pi_epsi Gamma

%----------------------------------------------------------------
% Steady State
%----------------------------------------------------------------

save SS C_SS Y_SS L_SS W_SS R_SS A_SS Z_SS B_SS Trans_SS MC_SS Pi_SS R_n_SS ...
    Wedge_EE_SS c_opt_SS ap_opt_SS mu_tilde_SS lambda_SS lambda_orig_SS

%----------------------------------------------------------------
% Other Quantities
%----------------------------------------------------------------

save aux grid_a grid_a_0 spliorder tol_GS_choice tol_GS_VFI tol_GS_dist ...
    states Phi_epsi Phi Emat fspaceerga fspace