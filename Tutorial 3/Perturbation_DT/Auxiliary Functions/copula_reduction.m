function lambda = copula_reduction(lambda_epsi,lambda_a,lambda_SS)

%----------------------------------------------------------------
% Size Inputs
%----------------------------------------------------------------

n_a_dist    = length(lambda_a);
n_epsi      = length(lambda_epsi);

%----------------------------------------------------------------
% Marginal CDFs
%----------------------------------------------------------------

% steady state

lambda_epsi_SS = sum(lambda_SS,2);
Lambda_epsi_SS = cumsum(lambda_epsi_SS); Lambda_epsi_SS(end) = 1;

lambda_a_SS    = sum(lambda_SS,1)';
Lambda_a_SS    = cumsum(lambda_a_SS); Lambda_a_SS(end) = 1;

Lambda_SS = cumsum(cumsum(lambda_SS,1),2); Lambda_SS(end,end) = 1;

% current values

Lambda_epsi = cumsum(lambda_epsi); Lambda_epsi(end) = 1;
Lambda_a    = cumsum(lambda_a); Lambda_a(end) = 1;

%----------------------------------------------------------------
% Get Joint CDF
%----------------------------------------------------------------

V = interpTwoD(Lambda_epsi,Lambda_a,Lambda_epsi_SS,Lambda_a_SS);
Lambda = reshape(V*Lambda_SS(:),n_epsi,n_a_dist);

%----------------------------------------------------------------
% Get Joint PDF
%----------------------------------------------------------------

lambda = NaN(n_epsi,n_a_dist);
for i_epsi = 1:n_epsi
    for i_a = 1:n_a_dist
        if i_epsi == 1 && i_a == 1
            lambda(i_epsi,i_a) = Lambda(i_epsi,i_a);
        elseif i_epsi > 1 && i_a == 1
            lambda(i_epsi,i_a) = Lambda(i_epsi,i_a) - Lambda(i_epsi-1,i_a);
        elseif i_epsi == 1 && i_a > 1
            lambda(i_epsi,i_a) = Lambda(i_epsi,i_a) - Lambda(i_epsi,i_a-1);
        else
            lambda(i_epsi,i_a) = Lambda(i_epsi,i_a) - Lambda(i_epsi-1,i_a) - Lambda(i_epsi,i_a-1) + Lambda(i_epsi-1,i_a-1);
        end
    end
end

end