pvec            = nodeunif(n_y_exp1,0.005,1-0.005);         % Make an equi-spaced grid in probabilities
e               = norminv(pvec,0,sigma_y);                % Invert normal for shocks
w               = normpdf(e,0,sigma_y);                   % Compute pdf of shocks
w               = w/sum(w);                                 % Normalise
fspaceZ         = fundef({'spli',grid_y_dist,0,1});              % Linear interpolant
gfun            = @(z,e) max(min(exp(rho_y*log(z)+e),max(grid_y_dist)),min(grid_y_dist));   % Constrained to lie within nodes
QZ              = zeros(n_s_dist,n_y_dist);
P               = zeros(n_y_dist,n_y_dist);                           % P constructed so can compute steady state Psszf and compare to Pssz
for i = 1:n_y_exp1
    g           = gfun(states_dist(:,2),e(i));
    QZi         = funbas(fspaceZ,g);
    QZ          = QZ + w(i)*QZi;
    P           = P  + w(i)*funbas(fspaceZ,gfun(grid_y_dist,e(i)));
end

clear pvec e w fspaceZ P g QZi