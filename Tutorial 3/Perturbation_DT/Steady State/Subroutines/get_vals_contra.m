%% UPPER BOUND

%----------------------------------------------------------------
% Set Distance
%----------------------------------------------------------------

dist_VFI = 1;

%----------------------------------------------------------------
% Complete Household Inputs
%----------------------------------------------------------------

R_k_SS   = b_guess(1);
R_SS     = R_k_SS + (1-delta);
KL_SS    = (R_k_SS/(alpha * Z_SS))^(-1/(1-alpha));
W_SS     = (1-alpha) * Z_SS * KL_SS^alpha;

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------  

% asset grid
                    
if grid_a_VFI_lin == 1
    if a_min < 0
        n_a_VFI_neg      = round(1/4 * n_a_VFI);
        n_a_VFI_pos      = n_a_VFI - n_a_VFI_neg;
        grid_a_VFI_neg_0 = linspace(a_min,0,n_a_VFI_neg);
        grid_a_VFI_pos_0 = linspace(0,a_max,n_a_VFI_pos);
        grid_a_VFI_0     = [grid_a_VFI_neg_0,grid_a_VFI_pos_0(2:end)];
    else
        grid_a_VFI_0 = linspace(a_min,a_max,n_a_VFI);
    end
else
    coeff_power = 0.9;
    power       = 8;
    if a_min < 0
        n_a_VFI_neg      = round(1/4 * n_a_VFI);
        n_a_VFI_pos      = n_a_VFI - n_a_VFI_neg;
        grid_a_VFI_neg_0 = linspace(a_min,0,n_a_VFI_neg);
        grid_a_VFI_pos_0 = linspace(0,1,n_a_VFI_pos);
        grid_a_VFI_pos_0 = 0 + (a_max-0)*((1 - coeff_power) * grid_a_VFI_pos_0 + coeff_power * (grid_a_VFI_pos_0.^power));
        grid_a_VFI_0     = [grid_a_VFI_neg_0,grid_a_VFI_pos_0(2:end)];
    else
        grid_a_VFI_0   	 = linspace(0,1,n_a_VFI_0);
        grid_a_VFI_0	 = a_min + (a_max-a_min)*((1 - coeff_power) * grid_a_VFI_0 + coeff_power * (grid_a_VFI_0.^power));
    end
end

%----------------------------------------------------------------
% Projection Approximation
%----------------------------------------------------------------

fspace          = fundef({'spli',grid_a_VFI_0,0,spliorder(1)},...
                         {'spli',grid_epsi,0,spliorder(2)});
states_VFI_grid  = funnode(fspace);
states_VFI       = gridmake(states_VFI_grid);

grid_a_VFI     = states_VFI(states_VFI(:,2)==states_VFI(1,2),1)';
n_a_VFI        = size(grid_a_VFI,2);

n_s_VFI  = n_a_VFI * n_epsi;

Phi_epsi_VFI   = splibas(grid_epsi,0,spliorder(2),states_VFI(:,2));
Phi_A_VFI      = splibas(grid_a_VFI_0,0,spliorder(1),states_VFI(:,1));
Phi_VFI        = dprod(Phi_epsi_VFI,Phi_A_VFI);
Emat_VFI       = kron(Pi_epsi,speye(n_a_VFI))*Phi_VFI;

%----------------------------------------------------------------
% Choice Regions
%----------------------------------------------------------------

optset('goldenx','tol',tol_GS_choice)

% lower bound of choice region for consumption

c_lb = 1e-8 * ones(n_s_VFI,1);

% upper bound of choice region for consumption

obj_ub  = @(c) -abs(c + a_min - (states_VFI(:,2) .* W_SS .* L_SS + R_SS * states_VFI(:,1)));
B_ub    = [zeros(n_s_VFI,1),2 * a_max * ones(n_s_VFI,1)];
c_ub    = goldenx(obj_ub,B_ub(:,1),B_ub(:,2));

%----------------------------------------------------------------
% VFI
%----------------------------------------------------------------

optset('goldenx','tol',tol_GS_VFI)

if R_it == 1
    h  = zeros(2*n_s_VFI,1);
else
    if use_prev == 0
       h  = zeros(2*n_s_VFI,1);
    end
end
h1 = h(1:n_s_VFI,1);
h2 = h(n_s_VFI+1:end,1);

while dist_VFI > tol_VFI
    obj_VFI = @(c) obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);
    c_opt   = goldenx(obj_VFI,c_lb,c_ub);
    for j = 1:totit_Howard
        v1       = obj_VFI(c_opt);
        h1_new   = Phi_VFI\v1; 
        h2_new   = Phi_VFI\(Emat_VFI * h1);
        h_new    = [h1_new;h2_new];
        h        = h_new;
        h1       = h(1:n_s_VFI,1);
        h2       = h(n_s_VFI+1:end,1);
        obj_VFI = @(c) obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);
    end
    obj_VFI = @(c) obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);
    if do_newton == 0
        c_opt    = goldenx(obj_VFI,c_lb,c_ub);
        v1       = obj_VFI(c_opt);
        h1_new   = Phi_VFI\v1; 
        h2_new   = Phi_VFI\(Emat_VFI * h1);
        h_new    = [h1_new;h2_new];
        dist_VFI = norm(h_new - h)/norm(h);
    elseif do_newton == 1
        c_opt      = goldenx(obj_VFI,c_lb,c_ub);
        ap_opt     = W_SS .* L_SS .* states_VFI(:,2) + R_SS * states_VFI(:,1) - c_opt;
        v1         = obj_VFI(c_opt);
        v2         = Emat_VFI * h1;
        Phi_Ap_opt = splibas(grid_a_VFI_0,0,spliorder(1),ap_opt);
        Phi_ApY    = dprod(Phi_epsi_VFI,Phi_Ap_opt);
        Jacobian_VFI = [ Phi_VFI,    -beta*Phi_ApY;
                        -Emat_VFI,              Phi_VFI];
        h_new        = h - Jacobian_VFI\([Phi_VFI*h1 - v1 ;
                                   Phi_VFI*h2 - v2]);
        dist_VFI = norm(h_new - h)/norm(h);
    end
    if disp_VFIdist == 1
       disp(dist_VFI)
    end
end

% collect results

h = h_new;

%----------------------------------------------------------------
% Distribution: Grids
%----------------------------------------------------------------

% asset grid

if grid_a_dist_lin == 1
    if a_min < 0
        n_a_dist_neg      = round(1/4 * n_a_dist);
        n_a_dist_pos      = n_a_dist - n_a_dist_neg;
        grid_a_dist_neg_0 = linspace(a_min,0,n_a_dist_neg);
        grid_a_dist_pos_0 = linspace(0,a_max,n_a_dist_pos);
        grid_a_dist_0     = [grid_a_dist_neg_0,grid_a_dist_pos_0(2:end)];
    else
        grid_a_dist_0 = linspace(a_min,a_max,n_a_dist);
    end
else
    coeff_power = 0.9;
    power       = 8;
    if a_min < 0
        n_a_dist_neg      = round(1/4 * n_a_dist);
        n_a_dist_pos      = n_a_dist - n_a_dist_neg;
        grid_a_dist_neg_0 = linspace(a_min,0,n_a_dist_neg);
        grid_a_dist_pos_0 = linspace(0,1,n_a_dist_pos);
        grid_a_dist_pos_0 = 0 + (a_max-0)*((1 - coeff_power) * grid_a_dist_pos_0 + coeff_power * (grid_a_dist_pos_0.^power));
        grid_a_dist_0     = [grid_a_dist_neg_0,grid_a_dist_pos_0(2:end)];
    else
        grid_a_dist_0   	 = linspace(0,1,n_a_dist_0);
        grid_a_dist_0	 = a_min + (a_max-a_min)*((1 - coeff_power) * grid_a_dist_0 + coeff_power * (grid_a_dist_0.^power));
    end
end

%----------------------------------------------------------------
% Distribution: Projection Approximation
%----------------------------------------------------------------

fspace          = fundef({'spli',grid_a_dist_0,0,spliorder(1)},...
                         {'spli',grid_epsi,0,spliorder(2)});
sgrid_dist      = funnode(fspace);
states_dist     = gridmake(sgrid_dist);

grid_a_dist      = states_dist((states_dist(:,2)==states_dist(1,2)),1)';
n_a_dist         = size(grid_a_dist,2);

n_s_dist  = n_a_dist * n_epsi;

Phi_epsi_dist   = splibas(grid_epsi,0,spliorder(2),states_dist(:,2));
Phi_A_dist      = splibas(grid_a_dist_0,0,spliorder(1),states_dist(:,1));
Phi_dist        = dprod(Phi_epsi_dist,Phi_A_dist);

%----------------------------------------------------------------
% Distribution: Optimal Choices
%----------------------------------------------------------------

% objective function

obj_dist = @(c) obj_dist_fun(c,states_dist,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_dist,h2);

% choice regions

optset('goldenx','tol',tol_GS_choice)

c_lb_dist = 1e-8 * ones(n_s_dist,1);

obj_ub_dist  = @(c) -abs(c + a_min - (states_dist(:,2) .* W_SS .* L_SS ...
                    + R_SS * states_dist(:,1)));
B_ub_dist    = [zeros(n_s_dist,1),2 * a_max * ones(n_s_dist,1)];
c_ub_dist    = goldenx(obj_ub_dist,B_ub_dist(:,1),B_ub_dist(:,2));
    
% optimal choices
   
c_opt_dist = goldenx(obj_dist,c_lb_dist,c_ub_dist);
ap_opt_dist = W_SS .* L_SS .* states_dist(:,2) + R_SS * states_dist(:,1) - c_opt_dist;

%----------------------------------------------------------------
% Distribution: Final Computation
%----------------------------------------------------------------

ap_opt_dist   = max(min(ap_opt_dist,a_max),a_min);
fspaceerga    = fundef({'spli',grid_a_dist,0,1});

QZ = kron(Pi_epsi,ones(n_a_dist,1));

QA = funbas(fspaceerga,ap_opt_dist);
Q  = dprod(QZ,QA);

lambda_SS      = full(ergodicdist(Q));
lambda_orig_SS = lambda_SS;
lambda_SS      = permute(reshape(lambda_SS,[n_a_dist,n_epsi]),[2 1]);

%----------------------------------------------------------------
% Compute Aggregates
%----------------------------------------------------------------

% get re-order policy function

c_opt_SS  = permute(reshape(c_opt_dist,[n_a_dist,n_epsi]),[2 1]);
ap_opt_SS = permute(reshape(ap_opt_dist,[n_a_dist,n_epsi]),[2 1]);

% compute other aggregates

A_grid        = NaN(n_epsi,n_a_dist); % asset holdings
for i_epsi = 1:n_epsi
    for i_a = 1:n_a_dist
            A_grid(i_epsi,i_a) = grid_a_dist(i_a);                   
    end
end

A_SS  = sum(sum(A_grid(:,:).*lambda_SS(:,:)));

%----------------------------------------------------------------
% Check Asset Market Clearing
%----------------------------------------------------------------

K_SS = KL_SS * L_SS;
Y_SS = Z_SS * K_SS^alpha * L_SS^(1-alpha);
A_err = A_SS - K_SS;

f_b_guess(1) = A_err;

%% LOWER BOUND

%----------------------------------------------------------------
% Set Distance
%----------------------------------------------------------------

dist_VFI = 1;

%----------------------------------------------------------------
% Complete Household Inputs
%----------------------------------------------------------------

R_k_SS   = a_guess(1);
R_SS     = R_k_SS + (1-delta);
KL_SS    = (R_k_SS/(alpha * Z_SS))^(-1/(1-alpha));
W_SS     = (1-alpha) * Z_SS * KL_SS^alpha;

%----------------------------------------------------------------
% Grids
%----------------------------------------------------------------  

% asset grid
                    
if grid_a_VFI_lin == 1
    if a_min < 0
        n_a_VFI_neg      = round(1/4 * n_a_VFI);
        n_a_VFI_pos      = n_a_VFI - n_a_VFI_neg;
        grid_a_VFI_neg_0 = linspace(a_min,0,n_a_VFI_neg);
        grid_a_VFI_pos_0 = linspace(0,a_max,n_a_VFI_pos);
        grid_a_VFI_0     = [grid_a_VFI_neg_0,grid_a_VFI_pos_0(2:end)];
    else
        grid_a_VFI_0 = linspace(a_min,a_max,n_a_VFI);
    end
else
    coeff_power = 0.9;
    power       = 8;
    if a_min < 0
        n_a_VFI_neg      = round(1/4 * n_a_VFI);
        n_a_VFI_pos      = n_a_VFI - n_a_VFI_neg;
        grid_a_VFI_neg_0 = linspace(a_min,0,n_a_VFI_neg);
        grid_a_VFI_pos_0 = linspace(0,1,n_a_VFI_pos);
        grid_a_VFI_pos_0 = 0 + (a_max-0)*((1 - coeff_power) * grid_a_VFI_pos_0 + coeff_power * (grid_a_VFI_pos_0.^power));
        grid_a_VFI_0     = [grid_a_VFI_neg_0,grid_a_VFI_pos_0(2:end)];
    else
        grid_a_VFI_0   	 = linspace(0,1,n_a_VFI_0);
        grid_a_VFI_0	 = a_min + (a_max-a_min)*((1 - coeff_power) * grid_a_VFI_0 + coeff_power * (grid_a_VFI_0.^power));
    end
end

%----------------------------------------------------------------
% Projection Approximation
%----------------------------------------------------------------

fspace          = fundef({'spli',grid_a_VFI_0,0,spliorder(1)},...
                         {'spli',grid_epsi,0,spliorder(2)});
states_VFI_grid  = funnode(fspace);
states_VFI       = gridmake(states_VFI_grid);

grid_a_VFI     = states_VFI(states_VFI(:,2)==states_VFI(1,2),1)';
n_a_VFI        = size(grid_a_VFI,2);

n_s_VFI  = n_a_VFI * n_epsi;

Phi_epsi_VFI   = splibas(grid_epsi,0,spliorder(2),states_VFI(:,2));
Phi_A_VFI      = splibas(grid_a_VFI_0,0,spliorder(1),states_VFI(:,1));
Phi_VFI        = dprod(Phi_epsi_VFI,Phi_A_VFI);
Emat_VFI       = kron(Pi_epsi,speye(n_a_VFI))*Phi_VFI;

%----------------------------------------------------------------
% Choice Regions
%----------------------------------------------------------------

optset('goldenx','tol',tol_GS_choice)

% lower bound of choice region for consumption

c_lb = 1e-8 * ones(n_s_VFI,1);

% upper bound of choice region for consumption

obj_ub  = @(c) -abs(c + a_min - (states_VFI(:,2) .* W_SS .* L_SS + R_SS * states_VFI(:,1)));
B_ub    = [zeros(n_s_VFI,1),2 * a_max * ones(n_s_VFI,1)];
c_ub    = goldenx(obj_ub,B_ub(:,1),B_ub(:,2));

%----------------------------------------------------------------
% VFI
%----------------------------------------------------------------

optset('goldenx','tol',tol_GS_VFI)

if R_it == 1
    h  = zeros(2*n_s_VFI,1);
else
    if use_prev == 0
       h  = zeros(2*n_s_VFI,1);
    end
end
h1 = h(1:n_s_VFI,1);
h2 = h(n_s_VFI+1:end,1);

while dist_VFI > tol_VFI
    obj_VFI = @(c) obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);
    c_opt   = goldenx(obj_VFI,c_lb,c_ub);
    for j = 1:totit_Howard
        v1       = obj_VFI(c_opt);
        h1_new   = Phi_VFI\v1; 
        h2_new   = Phi_VFI\(Emat_VFI * h1);
        h_new    = [h1_new;h2_new];
        h        = h_new;
        h1       = h(1:n_s_VFI,1);
        h2       = h(n_s_VFI+1:end,1);
        obj_VFI = @(c) obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);
    end
    obj_VFI = @(c) obj_VFI_fun(c,states_VFI,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_VFI,h2);
    if do_newton == 0
        c_opt    = goldenx(obj_VFI,c_lb,c_ub);
        v1       = obj_VFI(c_opt);
        h1_new   = Phi_VFI\v1; 
        h2_new   = Phi_VFI\(Emat_VFI * h1);
        h_new    = [h1_new;h2_new];
        dist_VFI = norm(h_new - h)/norm(h);
    elseif do_newton == 1
        c_opt      = goldenx(obj_VFI,c_lb,c_ub);
        ap_opt     = W_SS .* L_SS .* states_VFI(:,2) + R_SS * states_VFI(:,1) - c_opt;
        v1         = obj_VFI(c_opt);
        v2         = Emat_VFI * h1;
        Phi_Ap_opt = splibas(grid_a_VFI_0,0,spliorder(1),ap_opt);
        Phi_ApY    = dprod(Phi_epsi_VFI,Phi_Ap_opt);
        Jacobian_VFI = [ Phi_VFI,    -beta*Phi_ApY;
                        -Emat_VFI,              Phi_VFI];
        h_new        = h - Jacobian_VFI\([Phi_VFI*h1 - v1 ;
                                   Phi_VFI*h2 - v2]);
        dist_VFI = norm(h_new - h)/norm(h);
    end
    if disp_VFIdist == 1
       disp(dist_VFI)
    end
end

% collect results

h = h_new;

%----------------------------------------------------------------
% Distribution: Grids
%----------------------------------------------------------------

% asset grid

if grid_a_dist_lin == 1
    if a_min < 0
        n_a_dist_neg      = round(1/4 * n_a_dist);
        n_a_dist_pos      = n_a_dist - n_a_dist_neg;
        grid_a_dist_neg_0 = linspace(a_min,0,n_a_dist_neg);
        grid_a_dist_pos_0 = linspace(0,a_max,n_a_dist_pos);
        grid_a_dist_0     = [grid_a_dist_neg_0,grid_a_dist_pos_0(2:end)];
    else
        grid_a_dist_0 = linspace(a_min,a_max,n_a_dist);
    end
else
    coeff_power = 0.9;
    power       = 8;
    if a_min < 0
        n_a_dist_neg      = round(1/4 * n_a_dist);
        n_a_dist_pos      = n_a_dist - n_a_dist_neg;
        grid_a_dist_neg_0 = linspace(a_min,0,n_a_dist_neg);
        grid_a_dist_pos_0 = linspace(0,1,n_a_dist_pos);
        grid_a_dist_pos_0 = 0 + (a_max-0)*((1 - coeff_power) * grid_a_dist_pos_0 + coeff_power * (grid_a_dist_pos_0.^power));
        grid_a_dist_0     = [grid_a_dist_neg_0,grid_a_dist_pos_0(2:end)];
    else
        grid_a_dist_0   	 = linspace(0,1,n_a_dist_0);
        grid_a_dist_0	 = a_min + (a_max-a_min)*((1 - coeff_power) * grid_a_dist_0 + coeff_power * (grid_a_dist_0.^power));
    end
end

%----------------------------------------------------------------
% Distribution: Projection Approximation
%----------------------------------------------------------------

fspace          = fundef({'spli',grid_a_dist_0,0,spliorder(1)},...
                         {'spli',grid_epsi,0,spliorder(2)});
sgrid_dist      = funnode(fspace);
states_dist     = gridmake(sgrid_dist);

grid_a_dist      = states_dist((states_dist(:,2)==states_dist(1,2)),1)';
n_a_dist         = size(grid_a_dist,2);

n_s_dist  = n_a_dist * n_epsi;

Phi_epsi_dist   = splibas(grid_epsi,0,spliorder(2),states_dist(:,2));
Phi_A_dist      = splibas(grid_a_dist_0,0,spliorder(1),states_dist(:,1));
Phi_dist        = dprod(Phi_epsi_dist,Phi_A_dist);

%----------------------------------------------------------------
% Distribution: Optimal Choices
%----------------------------------------------------------------

% objective function

obj_dist = @(c) obj_dist_fun(c,states_dist,L_SS,R_SS,W_SS,gamma,beta,grid_a_VFI_0,spliorder,Phi_epsi_dist,h2);

% choice regions

optset('goldenx','tol',tol_GS_choice)

c_lb_dist = 1e-8 * ones(n_s_dist,1);

obj_ub_dist  = @(c) -abs(c + a_min - (states_dist(:,2) .* W_SS .* L_SS ...
                    + R_SS * states_dist(:,1)));
B_ub_dist    = [zeros(n_s_dist,1),2 * a_max * ones(n_s_dist,1)];
c_ub_dist    = goldenx(obj_ub_dist,B_ub_dist(:,1),B_ub_dist(:,2));
    
% optimal choices
   
c_opt_dist = goldenx(obj_dist,c_lb_dist,c_ub_dist);
ap_opt_dist = W_SS .* L_SS .* states_dist(:,2) + R_SS * states_dist(:,1) - c_opt_dist;

%----------------------------------------------------------------
% Distribution: Final Computation
%----------------------------------------------------------------

ap_opt_dist   = max(min(ap_opt_dist,a_max),a_min);
fspaceerga    = fundef({'spli',grid_a_dist,0,1});

QZ = kron(Pi_epsi,ones(n_a_dist,1));

QA = funbas(fspaceerga,ap_opt_dist);
Q  = dprod(QZ,QA);

lambda_SS      = full(ergodicdist(Q));
lambda_orig_SS = lambda_SS;
lambda_SS      = permute(reshape(lambda_SS,[n_a_dist,n_epsi]),[2 1]);

%----------------------------------------------------------------
% Compute Aggregates
%----------------------------------------------------------------

% get re-order policy function

c_opt_SS  = permute(reshape(c_opt_dist,[n_a_dist,n_epsi]),[2 1]);
ap_opt_SS = permute(reshape(ap_opt_dist,[n_a_dist,n_epsi]),[2 1]);

% compute other aggregates

A_grid        = NaN(n_epsi,n_a_dist); % asset holdings
for i_epsi = 1:n_epsi
    for i_a = 1:n_a_dist
            A_grid(i_epsi,i_a) = grid_a_dist(i_a);                   
    end
end

A_SS  = sum(sum(A_grid(:,:).*lambda_SS(:,:)));

%----------------------------------------------------------------
% Check Asset Market Clearing
%----------------------------------------------------------------

K_SS = KL_SS * L_SS;
Y_SS = Z_SS * K_SS^alpha * L_SS^(1-alpha);
A_err = A_SS - K_SS;

f_a_guess(1) = A_err;