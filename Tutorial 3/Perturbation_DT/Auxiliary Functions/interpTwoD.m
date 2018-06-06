function V = interpTwoD(x_fine,y_fine,x_knot,y_knot)
% Compute Linear interpolation over 2d grid. It extrapolates outside of
% knot points
%
% by SeHyoun Ahn, March 2017
%
% PARAMETERS:
%    x_fine = x grid to interpolate to
%    y_fine = y grid to interpolate to
%    x_knot = x knot points to interpolate from
%    y_knot = y knot points to interpolate from
%
% OUTPUT:
%    V = matrix giving interpolation values
%
% EXAMPLE:
%    x_fine = linspace(0,2,80)';
%    y_fine = linspace(-1,1,100)';
%    x_knot = linspace(0,1,10)';
%    y_knot = linspace(0,1,10)';
%
%    V = interpTwoD(x_fine,y_fine,x_knot,y_knot);
%
%    z = bsxfun(@plus,x_knot.^3+exp(-x_knot),(y_knot'-0.5).^2);
%    surf(x_fine,y_fine,reshape(V*z(:),80,100)');
%


loc_x = sum(bsxfun(@minus,x_fine,x_knot')>=0,2);
loc_x = min(loc_x,length(x_knot)-1);
loc_x = max(loc_x,1);

loc_y = sum(bsxfun(@minus,y_fine,y_knot')>=0,2);
loc_y = min(loc_y,length(y_knot)-1);
loc_y = max(loc_y,1);
loc_y = loc_y';

t_x = (x_fine-x_knot(loc_x))./(x_knot(loc_x+1)-x_knot(loc_x));
t_y = (y_fine-y_knot(loc_y))./(y_knot(loc_y+1)-y_knot(loc_y));
t_y = t_y';

ind_fine = 1:length(x_fine)*length(y_fine);
n_x = length(x_knot);

i = repmat(ind_fine',4,1);

j_SW = bsxfun(@plus,loc_x,n_x*(loc_y-1));
j_SE = bsxfun(@plus,loc_x+1,n_x*(loc_y-1));
j_NW = bsxfun(@plus,loc_x,n_x*(loc_y));
j_NE = bsxfun(@plus,loc_x+1,n_x*(loc_y));

v_SW = bsxfun(@times,(1-t_x),(1-t_y));
v_SE = bsxfun(@times,(t_x),(1-t_y));
v_NW = bsxfun(@times,(1-t_x),(t_y));
v_NE = bsxfun(@times,(t_x),(t_y));

j = [j_SW(:);j_SE(:);j_NW(:);j_NE(:)];
v = [v_SW(:);v_SE(:);v_NW(:);v_NE(:)];
V = sparse(i,j,v,length(x_fine)*length(y_fine),n_x*length(y_knot));