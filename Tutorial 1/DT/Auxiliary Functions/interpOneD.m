function V = interpOneD(x,x_knot)
% Compute Linear interpolation over 1d grid. It extrapolates outside of
% knot points
%
% by SeHyoun Ahn, March 2017
%
% PARAMETERS:
%    x_fine = x grid to interpolate to
%    x_knot = x knot points to interpolate from
%
% OUTPUT:
%    V = matrix giving interpolation values
%
% EXAMPLE:
%    x_fine = linspace(-1,1,300)';
%    x_knot = linspace(0,2,20)';
%    V = interpOneD(x_fine,x_knot);
%    plot(V*(x_knot.^5+x_knot*0.5));
%

% loc = sum(bsxfun(@minus,x,x_knot')>=0,2);
n_x = size(x,1);
loc = sum((repmat(x,1,n_x) - repmat(x_knot',n_x,1))>=0,2);
loc = min(loc,length(x_knot)-1);
loc = max(loc,1);

t = (x-x_knot(loc))./(x_knot(loc+1)-x_knot(loc));
ind_x = 1:length(x);

i = repmat(ind_x',2,1);
j = [loc;loc+1];
v = [(1-t);t];

V = sparse(i,j,v,length(x),length(x_knot));