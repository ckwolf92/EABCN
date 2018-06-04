function quantile = quantile_fun(grid,dist,quantile_val)

lb_indx = find(cumsum(dist) > quantile_val,1);
ub_indx = lb_indx + 1;
dist_cumsum = cumsum(dist);
lb_prob = dist_cumsum(lb_indx);
ub_prob = dist_cumsum(ub_indx);
lb_weight = (ub_prob - quantile_val)/(ub_prob - lb_prob);
lb_val = grid(lb_indx);
ub_val = grid(ub_indx);
quantile = lb_weight * lb_val + (1-lb_weight) * ub_val;

end