function top_share = top_share_fun(grid,dist,quantile);

dist_cumsum = cumsum(dist,'reverse');
cutoff = max(find(dist_cumsum > quantile)+1);
weights = 0 * grid;
weights(cutoff:end) = 1;
weights(cutoff-1) = (quantile - dist_cumsum(cutoff))/(dist_cumsum(cutoff-1) - dist_cumsum(cutoff));
top_share = sum(grid .* dist .* weights)/sum(grid .* dist);

end