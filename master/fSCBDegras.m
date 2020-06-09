function z = fSCBDegras(covariance_onb,onbMatrix)

[nGridSpace, nBasis] = size( onbMatrix );

% make sure it's symetric (potential numerical nonsymetry)
covariance_onb = real(covariance_onb+covariance_onb')/2;

% make sure it's pos definite
[V,D] = eig(covariance_onb);
d = diag(D);

% if there's a low eigenvalue, write a warning
if (min(D(:)) < -0.1)
    fprintf('Warning: fSCBDegras: matrix is far away from being pos-def\n')
    d
end

d = max(d,1e-5); % make is pos-definite
D = diag(d); % back to a diagonal matrix
covariance_onb = V * D * inv(V); % make it positive definite

% the following is to calculate the correlation matrix in ONB
covariance_space = onbMatrix * covariance_onb * onbMatrix';
correlation_space = covariance_space ./ ( sqrt(diag(covariance_space)) * sqrt(diag(covariance_space))' );
correlation_onb = onbMatrix' * correlation_space * onbMatrix / nGridSpace^2;

mcrl = mvnrnd(zeros(1000,nBasis), correlation_onb); % sample from the correlation
mcrl_max =  max(abs(onbMatrix*mcrl')); % for each trajectory calculate the max abs value
z = quantile(mcrl_max, 0.95); % the quantile

end
