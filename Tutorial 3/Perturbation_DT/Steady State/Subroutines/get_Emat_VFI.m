pvec            = nodeunif(n_y_exp1,0.005,1-0.005);     % Make an equi-spaced grid in probabilities
e               = norminv(pvec,0,sigma_y);            
w               = normpdf(e,0,sigma_y);               % Invert normal for shocks
w               = w/sum(w);                             % Compute pdf of shocks
iNe             = ones(n_y_exp1,1);                           
iNs             = ones(n_s_VFI,1);
gfun            = @(z,e) max(min(exp(rho_y*log(z)+e),max(grid_y_VFI)),min(grid_y_VFI));   % Constrained to lie within nodes
g               = gfun(kron(states_VFI(:,2),iNe),kron(iNs,e));
Phi             = funbas(fspace,[kron(states_VFI(:,1),iNe),g]);
Ikronw          = kron(eye(n_s_VFI),w');
Emat_VFI        = Ikronw*Phi;

clear pvec e w iNe iNs gfun g Phi Ikronw