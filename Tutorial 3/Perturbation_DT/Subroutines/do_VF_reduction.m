a_min = min(grid_a_VFI);
a_max = max(grid_a_VFI);

coeff_power = 0.9;
power       = 8;

knots   	 = linspace(0,1,n_VF_reduced);
knots	 = a_min + (a_max-a_min)*((1 - coeff_power) * knots + coeff_power * (knots.^power));

n_rest = n_g + n_x;

n_post = n_epsi;
n_prior = 1;
[from_spline_small, to_spline_small] = createQuadSplineNDa(grid_a_VFI',knots',n_prior,n_post,n_rest);

n_splined = size(from_spline_small,2)*(2*n_epsi);

to_spline_big = kron(speye(2*n_epsi),to_spline_small);
to_spline = spdiags(ones(n_splined+n_rest,1),n_v-n_splined,n_rest+n_splined,n_v+n_rest);
to_spline(1:n_splined,1:n_v) = to_spline_big;

from_spline_big = kron(speye(2*n_epsi),from_spline_small);
from_spline = spdiags(ones(n_splined+n_rest,1),n_splined-n_v,n_v+n_rest,n_rest+n_splined);
from_spline(1:n_v,1:n_splined) = from_spline_big;