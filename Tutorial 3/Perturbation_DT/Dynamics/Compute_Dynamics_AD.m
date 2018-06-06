%% HANK: COMPUTE SHOCK IRFs
% Christian Wolf
% this version: 03/14/2018

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/ckwolf/Dropbox/Research/Miscellaneous/EABCN_Chile/Slides_CKW/codes/Tutorial 3/Perturbation_DT';  
% addpath('/Users/ckwolf/Dropbox/Coding/compecon2011/CEtools')
% addpath('/Users/ckwolf/Dropbox/Coding/compecon2011/CEdemos')
  
addpath([path '/Auxiliary Functions'])
addpath([path '/Subroutines'])
addpath([path '/Steady State'])

cd([path '/Dynamics/']);

addpath('/Applications/Dynare/4.5.1/matlab/')

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

%----------------------------------------------------------------
% Imports
%----------------------------------------------------------------

load param_agg
load param_households
load SS
load aux

%----------------------------------------------------------------
% Some Auxiliary Quantities
%----------------------------------------------------------------

n_a_VFI  = size(grid_a,2);
n_a_dist = size(grid_a,2);
n_epsi   = size(grid_epsi,2);

n_s_VFI  = n_a_VFI * n_epsi;
n_s_dist = n_a_dist * n_epsi;

a_min = min(grid_a);
a_max = max(grid_a);

%----------------------------------------------------------------
% Solution Settings
%----------------------------------------------------------------

global reduce_VF reduce_dist

% value function

reduce_VF = 0;
VF_knots  = 5;

% distribution

reduce_dist = 0;

%% REDUCTION

%----------------------------------------------------------------
% Value Function
%----------------------------------------------------------------

global from_spline to_spline

if reduce_VF == 1
    get_VF_splines
else
    to_spline   = speye(n_s_VFI);
    from_spline = speye(n_s_VFI);
end

n_a_VFI_reduced = size(from_spline,2)/n_epsi;
n_s_VFI_reduced = size(from_spline,2);

%----------------------------------------------------------------
% Distribution
%----------------------------------------------------------------

global lambda_epsi_SS lambda_a_SS

if reduce_dist == 1
    n_s_dist_reduced = n_a_dist;
else
    n_s_dist_reduced = n_s_dist;
end

lambda_epsi_SS = sum(lambda_SS,2);
lambda_a_SS    = sum(lambda_SS,1)';

%% SET UP SYSTEM

%----------------------------------------------------------------
% System Size
%----------------------------------------------------------------

global n_shocks nVars nEErrors

n_v      = n_s_VFI_reduced;
n_g      = n_s_dist_reduced;
n_x      = 10;

n_shocks  = 1;
nVars     = n_v + n_g + n_x;
nEErrors  = n_s_VFI_reduced + 1;

%----------------------------------------------------------------
% Collect Steady State Objects
%----------------------------------------------------------------

global vars_SS

vars_SS     = NaN(nVars,1);

vars_SS(1:n_s_VFI) = mu_tilde_SS;

if reduce_dist == 1
    vars_SS(n_s_VFI+1:n_s_VFI+n_g) = lambda_a_SS;
else
    vars_SS(n_s_VFI+1:n_s_VFI+n_g) = lambda_orig_SS;
end

vars_SS(n_s_VFI+n_s_dist_reduced+1)  = C_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+2)  = Y_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+3)  = L_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+4)  = W_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+5)  = A_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+6)  = MC_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+7)  = Pi_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+8)  = R_n_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+9)  = Trans_SS;
vars_SS(n_s_VFI+n_s_dist_reduced+10) = 0;

%% LINEARIZATION

t0 = tic;
fprintf('Taking derivatives of equilibrium conditions...\n')

vars 					= zeros(nVars + nVars + nEErrors + n_shocks,1);
vars 					= myAD(vars);                               % converts to dual numbers, required for differentiation
derivativesIntermediate = equilibriumConditions_reduced_AD(vars);              % evaluates function and derivatives
derivs 					= getderivs(derivativesIntermediate);       % extracts only derivatives

% Translate to canonical LRE form
g1                      = -derivs(:,1:nVars);
g0                      = derivs(:,nVars+1:2*nVars);
pi                      = -derivs(:,2*nVars+1:2*nVars+nEErrors);
psi                     = -derivs(:,2*nVars+nEErrors+1:2*nVars+nEErrors+n_shocks);
constant				= sparse(nVars,1);

fprintf('Time to take derivatives: %2.4f seconds\n\n\n',toc(t0))

%% SOLUTION

t0 = tic;
fprintf('Solving linear system...\n')

[G1,~,impact,~,~,~,~,eu] = gensys(full(g0),full(g1),full(constant),full(psi),full(pi));

fprintf('...Done!\n')
fprintf('Existence and uniqueness? %2.0f and %2.0f\n',eu);
fprintf('Time to solve linear system: %2.4f seconds\n\n\n',toc(t0))

%% SIMULATION: IMPULSE RESPONSE

t0				= tic;
fprintf('Simulating IRFs...\n')

%----------------------------------------------------------------
% Preparations
%----------------------------------------------------------------

% Set parameters of the simulation
T                   = 100;
N_sim				= T;

% Preallocate matrices
vTime		= linspace(1,T,N_sim)';
[nVars,~]   = size(G1);
mVars       = zeros(nVars,N_sim,n_shocks);

%----------------------------------------------------------------
% Simulation
%----------------------------------------------------------------

for shock = 1:n_shocks    

% Vector of aggregate shocks
vAggregateShock				= zeros(N_sim,1);
vAggregateShock(1)          = 1;

% Simulate
for n = 1 : N_sim
	mVars(:,n+1,shock)	= G1 * mVars(:,n,shock) + impact(:,shock) * vAggregateShock(n,1)';	
end

end

%----------------------------------------------------------------
% Collect Results
%----------------------------------------------------------------

% remove crap

mVars     = mVars(:,2:end,:);
mVars_aux = mVars(1:n_v+n_g,:,:);
mVars     = mVars(n_v+n_g+1:end,:,:);

% monetary policy shock

c_e_m     = mVars(1,:,1);
y_e_m     = mVars(2,:,1);
l_e_m     = mVars(3,:,1);
w_e_m     = mVars(4,:,1);
a_e_m     = mVars(5,:,1);
mc_e_m    = mVars(6,:,1);
pi_e_m    = mVars(7,:,1);
r_n_e_m   = mVars(8,:,1);
trans_e_m = mVars(9,:,1);
u_m_e_m   = mVars(10,:,1);

fprintf('...Done!\n')
fprintf('Time to simulate model: %2.4f seconds\n\n\n',toc(t0))

%----------------------------------------------------------------
% Get Wedges
%----------------------------------------------------------------

wedge_EE_e_m = c_e_m(1:end-1) - (c_e_m(2:end) - 1/gamma * (r_n_e_m(1:end-1) - pi_e_m(2:end)));

%% PLOTS

T_plot = 50;
vTime  = 0:1:T_plot;

plot_results