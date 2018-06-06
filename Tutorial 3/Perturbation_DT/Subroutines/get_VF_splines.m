coeff_power = 0.9;
power       = 8;

knots = linspace(0,1,VF_knots);
knots = a_min + (a_max-a_min)*((1 - coeff_power) * knots + coeff_power * (knots.^power));

n_post = n_epsi;
n_prior = 1;
[from_spline_small, to_spline_small] = createQuadSplineNDa(grid_a',knots');

to_spline = kron(speye(n_epsi),to_spline_small);

from_spline = kron(speye(n_epsi),from_spline_small);