function [Pi_epsi, epsi] = rouwen_star(rho, sigma, n_epsi, p_star) 
    [Pi_epsi, epsi] = rouwen(rho,0,sigma,n_epsi);
    epsi            = epsi';
    Pi_epsi         = Pi_epsi';                         % row i, column j gives the probability of moving from epsi(i) to epsi(j)

%     epsi(5) = epsi(5)*(1+2*p_star) ; 
%     Pi_epsi(:,5) = Pi_epsi(:,5) + p_star ; 
%     Pi_epsi(1,2) = Pi_epsi(1,2) - p_star ; 
%     Pi_epsi(2,3) = Pi_epsi(2,3) - p_star ; 
%     Pi_epsi(3,4) = Pi_epsi(3,4) - 3/4*p_star ; 
%     Pi_epsi(3,2) = Pi_epsi(3,2) - 1/4*p_star ; 
%     Pi_epsi(4,3) = Pi_epsi(4,3) - p_star ; 
%     Pi_epsi(5,4) = Pi_epsi(5,4) - p_star ; 
    
end