%% HETEROGENEOUS-HOUSEHOLD MODELS: SOLVE INCOME-FLUCTUATION PROBLEM
% Christian Wolf
% this version: 05/26/2018

%% HOUSEKEEPING
 
clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/ckwolf/Dropbox/Research/Miscellaneous/EABCN_Chile/Slides_CKW/codes/Tutorial 1/DT';

addpath([path '/Auxiliary Functions'])
addpath([path '/Income Process'])
cd(path);

%% ECONOMIC PARAMETERS

%----------------------------------------------------------------
% Household
%----------------------------------------------------------------

% preferences

beta          = 0.98;
gamma         = 1;
probdeath     = 1/(4*40);
beta_hat      = beta * (1-probdeath);

% initial wealth

wealth_newborn = 0;

% income risk

load logyPgrid.txt
load yPtrans.txt
grid_yP = exp(logyPgrid);
n_yP = length(grid_yP);
Pi_yP = yPtrans;
clear logyPgrid yPtrans

load logyTgrid.txt
load yTdist.txt
grid_yT = exp(logyTgrid);
n_yT = length(grid_yT);
Pi_yT = repmat(yTdist',n_yT,1);
clear logyTgrid yTdist

n_y = n_yT * n_yP;
grid_y = repmat(grid_yP,1,n_yT) .* reshape(repmat(grid_yT,1,n_yP)',1,n_y);
indx_yP = reshape(repmat((1:1:n_yP)',1,n_yT),n_y,1)';
indx_yT = reshape(repmat(1:1:n_yT,n_yP,1),n_y,1)';

Pi_y = NaN(n_y,n_y);
for i_y = 1:n_y
    for i_yy = 1:n_y
        Pi_y(i_y,i_yy) = Pi_yP(indx_yP(i_y),indx_yP(i_yy)) * Pi_yT(indx_yT(i_y),indx_yT(i_yy));
    end
end

y_dist = ergodicdist(Pi_y)'; 
grid_y = grid_y/sum(grid_y(:).*y_dist(:));
mean_y = grid_y * y_dist';

%----------------------------------------------------------------
% Aggregates
%----------------------------------------------------------------

% interest rate

r_SS     = 0.01;
r_hat_SS = r_SS + 1/(1-probdeath) - 1;

% transfer for MPC distribution

gift = 0.01 * mean_y;

%% SOLUTION PARAMETERS

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------

n_a_0      = 50; % number of grid points 
a_min      = 0;
a_max      = 50 * max(grid_yP);
spliorder  = [3 1];

%----------------------------------------------------------------
% EGP Iteration
%----------------------------------------------------------------

EGP_tol      = 10^(-8);
disp_EGPdist = 1;

%% ASSEMBLE GRID

%----------------------------------------------------------------
% Raw Asset Grid
%----------------------------------------------------------------

coeff_power = 0.9;
power       = 8;
grid_a_0   	= linspace(0,1,n_a_0);
grid_a_0	= a_min + (a_max-a_min)*((1 - coeff_power) * grid_a_0 + coeff_power * (grid_a_0.^power));

%----------------------------------------------------------------
% Splines
%----------------------------------------------------------------

% put productivity and asset grids together

n_s_0  = n_a_0 * n_y;

fspace          = fundef({'spli',grid_a_0,0,spliorder(1)},...
                         {'spli',grid_y,0,spliorder(2)});
states_grid  = funnode(fspace);
states       = gridmake(states_grid);

grid_a        = states(states(:,2)==states(1,2),1)';
n_a           = size(grid_a,2);
[~,newborn_pos] = min(abs(grid_a - wealth_newborn));

n_s  = n_a * n_y;

Phi_epsi   = splibas(grid_y,0,spliorder(2),states(:,2));
Phi_A      = splibas(grid_a_0,0,spliorder(1),states(:,1));
Phi        = dprod(Phi_epsi,Phi_A);
Emat       = kron(Pi_y,speye(n_a));
Emat_yP    = Emat(1:n_a*n_yP,:);

%% MAIN SOLUTION LOOP  
 
%----------------------------------------------------------------
% Endogenous Gridpoint Iteration
%----------------------------------------------------------------

% preparations

dist_EGP = 1;
EGP_it   = 0;

% initial guess

cp_opt  = states(:,2) + r_hat_SS * states(:,1);
mup_opt = 0 * cp_opt;
mutilde_opt = Emat_yP * (cp_opt.^(-gamma));

% iteration

while dist_EGP > EGP_tol
    
% one step for EGP

[c_opt,ap_opt,mutilde_upd] = EGP_fun(mutilde_opt,r_hat_SS,1,1,0,beta_hat,gamma,0,states,grid_a,n_yT,Emat_yP);
    
% update

dist_EGP = norm(mutilde_upd - mutilde_opt)/norm(mutilde_opt);

if disp_EGPdist == 1 && mod(EGP_it,100) == 0
    disp(dist_EGP)
end

mutilde_opt = mutilde_upd;
EGP_it      = EGP_it + 1;

end

%----------------------------------------------------------------
% Distribution
%----------------------------------------------------------------

ap_opt     = max(min(ap_opt,a_max),a_min);
fspaceerga = fundef({'spli',grid_a,0,1});

QZ_live = kron(Pi_y,ones(n_a,1));
QA_live = funbas(fspaceerga,ap_opt);
Q_live  = dprod(QZ_live,QA_live);

QA_death = sparse(n_s,n_a);
QA_death(:,newborn_pos) = 1;
QZ_death = repmat(y_dist,n_s,1);
Q_death  = dprod(QZ_death,QA_death);

Q = (1-probdeath) * Q_live + probdeath * Q_death;

lambda_SS      = full(ergodicdist(Q));
lambda_orig_SS = lambda_SS;
lambda_SS      = permute(reshape(lambda_SS,[n_a,n_y]),[2 1]);

%----------------------------------------------------------------
% Compute Aggregates
%----------------------------------------------------------------

% re-order policy function

c_opt_SS  = permute(reshape(c_opt,[n_a,n_y]),[2 1]);
ap_opt_SS = permute(reshape(ap_opt,[n_a,n_y]),[2 1]);

% compute other aggregates

C_grid = c_opt_SS;
C_SS   = sum(sum(C_grid(:,:).*lambda_SS(:,:)));

A_grid = repmat(grid_a,n_y,1);
A_SS   = sum(sum(A_grid(:,:).*lambda_SS(:,:)));

%% STEADY-STATE STATISTICS: COMPUTATIONS

%----------------------------------------------------------------
% MPC
%----------------------------------------------------------------

[c_opt_MPC,~,~] = EGP_fun(mutilde_opt,r_hat_SS,1,1,0+gift,beta_hat,gamma,0,states,grid_a,n_yT,Emat_yP);

c_opt_MPC = permute(reshape(c_opt_MPC,[n_a,n_y]),[2 1]);
MPC_dist_SS = (c_opt_MPC - c_opt_SS)./gift;
MPC_avg_SS = sum(MPC_dist_SS(:) .* lambda_SS(:));

%% STEADY-STATE STATISTICS: PLOTS

%----------------------------------------------------------------
% MPC Distribution
%----------------------------------------------------------------

if min(MPC_dist_SS(:)) < -0.01
    warning('There are negative MPCs!')
end

MPC_grid = (0:0.05:1);
MPC_grid_prob = NaN(1,length(MPC_grid)-1);
for i_MPC = 1:length(MPC_grid)-1
    MPC_LB = MPC_grid(i_MPC);
    MPC_UB = MPC_grid(i_MPC+1);
    if i_MPC == length(MPC_grid)-1
        indic = (MPC_dist_SS >= MPC_LB);
    else
        indic = (MPC_dist_SS >= MPC_LB) & (MPC_dist_SS < MPC_UB);
    end
    MPC_grid_prob(i_MPC) = sum(indic(:) .* lambda_SS(:));
end

figure(1)
set(gca,'FontSize',16)
bar(MPC_grid(2:end),MPC_grid_prob)
xlabel('MPC')
ylabel('Mass')
xlim([-1 1])
title('MPC Distribution')

%----------------------------------------------------------------
% Policy Function
%----------------------------------------------------------------

y_plot = [1 6 13 17 21 28 33];
a_plot = round(n_a/2);

figure(4)
set(gca,'FontSize',16)
plot(grid_a(1:a_plot),squeeze(c_opt_SS(y_plot,(1:a_plot))))
xlabel('Assets, a')
title('Consumption Policy Functions')