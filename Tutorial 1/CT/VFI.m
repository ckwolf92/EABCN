%% HETEROGENEOUS-HOUSEHOLD MODELS: SOLVE INCOME-FLUCTUATION PROBLEM
% Christian Wolf
% this version: 05/26/2018

%% HOUSEKEEPING

clc
clear all
close all

warning('off','MATLAB:dispatcher:nameConflict')

path = '/Users/ckwolf/Dropbox/Research/Miscellaneous/EABCN_Chile/Slides_CKW/codes/Tutorial 1/CT';

addpath([path '/Income Process'])
cd(path);

%% ECONOMIC PARAMETERS

%----------------------------------------------------------------
% Economic Parameters 
%----------------------------------------------------------------

% preferences

rho       = 0.02;
gamma     = 2;
probdeath = 0;
rho_hat   = rho + probdeath;

% income risk

load ygrid_combined.txt
load ymarkov_combined.txt

grid_y = exp(ygrid_combined);
n_y = length(grid_y);

AT = ymarkov_combined';
dist_y = ones(n_y,1) / n_y;
for i_HJB = 1 : 50
	g_y_new = (speye(n_y) - AT * 1000) \ dist_y;
	diff = max(abs(g_y_new - dist_y));
	if diff < 10e-6
		break
	end
	dist_y = g_y_new;
end
grid_y = grid_y/sum(grid_y(:).*dist_y(:));
mean_y = sum(grid_y' * dist_y);
grid_y = grid_y';

%----------------------------------------------------------------
% Aggregates
%----------------------------------------------------------------

% interest rate

r_SS     = 0.01;
r_hat_SS = probdeath + r_SS;

%% SOLUTION PARAMETERS

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------

% asset grid

n_a   = 100;
a_min = 0;
a_max = 100 * mean_y;

coeff_power = 0.9; power = 8;
grid_a = linspace(0,1,n_a);
grid_a = a_min + (a_max-a_min)*((1 - coeff_power) * grid_a + coeff_power * (grid_a.^power));
grid_a = grid_a';

% auxiliary objects

Aswitch = kron(ymarkov_combined,speye(n_a));
grid_y_rep = ones(n_a,1) * grid_y;
grid_a_rep = grid_a * ones(1,n_y);

%----------------------------------------------------------------
% Solution Numerics
%----------------------------------------------------------------

maxit_HJB = 100;
crit_HJB = 1e-6;
Delta = 1000;

%% MAIN SOLUTION LOOP

%----------------------------------------------------------------
% Housekeeping
%----------------------------------------------------------------

% Initialize matrices for finite differences

dVbf = zeros(n_a,n_y);
dVbb = zeros(n_a,n_y);
dVzf = zeros(n_a,n_y);
dVzb = zeros(n_a,n_y);
dVzz = zeros(n_a,n_y);
c    = zeros(n_a,n_y);

% Create grids for finite differences scheme

daf = ones(n_a,1);
dab = ones(n_a,1);
daf(1:n_a-1) = grid_a(2:n_a) - grid_a(1:n_a-1);
dab(2:n_a) = grid_a(2:n_a) - grid_a(1:n_a-1);
daf(n_a) = daf(n_a-1); dab(1) = dab(2);
daa_f = daf * ones(1,n_y);
daa_b = dab * ones(1,n_y);

%----------------------------------------------------------------
% HJB Iteration
%----------------------------------------------------------------

% Initial guess of value function

if gamma == 1
    v = 1/rho_hat * log(grid_y_rep + r_hat_SS .* grid_a_rep);
else    
    v = ((grid_y_rep + r_hat_SS .* grid_a_rep) .^ (1 - gamma)) / (rho_hat * (1 - gamma));
end

% iteraton

for i_HJB = 1:maxit_HJB

    V = v;
    V_n(:,:,i_HJB) = V;

    % Compute forward difference
    dVf(1:n_a-1,:) = (V(2:n_a,:) - V(1:n_a-1,:)) ./ (grid_a_rep(2:n_a,:) - grid_a_rep(1:n_a-1,:));
    dVf(n_a,:) = (grid_y + r_hat_SS .* a_max) .^ (-gamma);	% will never be used, but impose a <= amax just in case

    % Compute backward difference
    dVb(2:n_a,:) = (V(2:n_a,:) - V(1:n_a-1,:)) ./ (grid_a_rep(2:n_a,:) - grid_a_rep(1:n_a-1,:));
    dVb(1,:) = (grid_y + r_hat_SS .* a_min) .^ (-gamma);	% state constraint boundary condition again

    % Consumption and savings with forward difference
    cf = max(dVf,10^(-8)) .^ (-1 / gamma);
    cf(n_a,:) = min(cf(n_a,:),grid_y_rep(n_a,:) + r_hat_SS.*grid_a_rep(n_a,:));
    sf = grid_y_rep + r_hat_SS.*grid_a_rep - cf;

    % Consumption and savings with backward difference
    cb = max(dVb,10^(-8)) .^ (-1 / gamma);
    sb = grid_y_rep + r_hat_SS.*grid_a_rep - cb;
    sb(1,:) = 0;

    % Consumption with no drift
    c0 = grid_y_rep + r_hat_SS.*grid_a_rep;
    dV0 = (c0) .^ (-gamma);

    % dV_Upwind makes choice between forward or backward difference based on sign of drift
    If = sf > 0;	% positive drift --> forward difference
    Ib = (sb < 0) .* (1 - If);	% negative drift --> backward difference
    I0 = 1 - If - Ib; % no drift

    % Consumption choice with upwind differences
    c = If .* cf + Ib .* cb + I0 .* c0;
    if gamma == 1
        u = log(c);
    else
        u = ((c) .^ (1 - gamma)) ./ (1 - gamma);
    end
    dV_Upwind = dVf .* If + dVb .* Ib + dV0 .* I0;
    savingsSS = grid_y_rep + r_hat_SS.*grid_a_rep - c;

    % Construct matrix for implicit updating scheme
    X = -min(sb,0) ./ daa_b;
    Y = -max(sf,0) ./ daa_f + min(sb,0) ./ daa_b;
    Z = max(sf,0) ./ daa_f;

    % The following is needed because of a peculiarity of spdiags
    updiag = [0];
    for j=1:n_y
        updiag=[updiag;Z(1:n_a-1,j);0];
    end
    centdiag = reshape(Y,n_a*n_y,1);
    lowdiag=X(2:n_a,1);
    for j=2:n_y
        lowdiag=[lowdiag;0;X(2:n_a,j)];
    end

    % Create A matrix
    AA=spdiags(centdiag,0,n_a*n_y,n_a*n_y)+spdiags([updiag;0],1,n_a*n_y,n_a*n_y)+spdiags([lowdiag;0],-1,n_a*n_y,n_a*n_y);
    A = AA + Aswitch;
    if max(abs(sum(A,2)))>10^(-8)
       disp('Improper Transition Matrix')
       break
    end

    % Solve for new value function
    B = (1/Delta + rho_hat)*speye(n_a*n_y) - A;

    u_stacked = reshape(u,n_a*n_y,1);
    V_stacked = reshape(V,n_a*n_y,1);
    vec = u_stacked + V_stacked/Delta;

    V_stacked = B\vec; % Implicit scheme

    V = reshape(V_stacked,n_a,n_y);
    Vchange = V - v;
    v = V; 

    % Update convergence criterion
    dist(i_HJB) = max(max(abs(Vchange)));
    if dist(i_HJB)<crit_HJB
        %disp('Value Function Converged, Iteration = ')
        %disp(n)
        break

    end

end

%----------------------------------------------------------------
% Stationary Distribution
%----------------------------------------------------------------

% Adjust to ensure mass preservation with non-uniform grid
da_tilde = 0.5*(dab + daf);
da_tilde(1) = daf(1); da_tilde(n_a) = dab(n_a);
da_stacked = reshape(da_tilde*ones(1,n_y),n_a*n_y,1);
grid_diag = spdiags(da_stacked,0,n_a*n_y,n_a*n_y);
A_adj = grid_diag*A*grid_diag^(-1); %mass preservation

% Normalization so pdf integrates to 1 
llambda0 = ones(n_a*n_y,1);
lambda_sum = llambda0'*ones(n_a*n_y,1)*da_stacked;
llambda0 = llambda0./lambda_sum;

% Adjust for death
aux = sparse(n_a,n_a);
aux(:,1)=probdeath*da_tilde/da_tilde(1);
aux2 = kron(speye(n_y),aux);
A_adj = A_adj + aux2 - probdeath*speye(n_a*n_y);

% Solve linear system for distribution using implicit method
AT = A_adj';
llambda = llambda0;
for i_HJB = 1:50
    gg_new = (speye(n_a*n_y) - AT*Delta)\llambda;
    g_sum_t(i_HJB) = llambda'*da_stacked;
    dist(i_HJB) = max(abs(gg_new-llambda));
    if dist(i_HJB)<10^(-6)
        break
    end
    llambda = gg_new;
end

% Store distribution
lambda_st = llambda;
A_st = (reshape(grid_a_rep,n_a*n_y,1).*lambda_st)'*da_stacked;
C = (reshape(c,n_a*n_y,1).*lambda_st)'*da_stacked;
c_st = c;
adot_st = grid_y_rep + r_hat_SS.*grid_a_rep - c;
clear gg

% Marginal distributions
lambda = reshape(lambda_st,n_a,n_y);
lambda_a = sum(lambda,2);
dist_y = da_tilde'*lambda;

%% STEADY-STATE STATISTICS: PLOTS

%----------------------------------------------------------------
% Policy Function
%----------------------------------------------------------------

y_plot = [1 6 13 17 21 28 33];
a_plot = round(n_a/2);

figure(1)
set(gca,'FontSize',16)
plot(grid_a(1:a_plot),squeeze(c((1:a_plot),y_plot)))
xlabel('Assets, a')
title('Consumption Policy Functions')

%----------------------------------------------------------------
% Policy Function
%----------------------------------------------------------------

y_plot = [1];
a_plot = n_a;

figure(2)
set(gca,'FontSize',16)
plot(grid_a(1:a_plot),squeeze(lambda((1:a_plot),y_plot)))
xlabel('Assets, a')
title('Wealth Distribution')