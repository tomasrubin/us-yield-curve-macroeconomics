function mu_est = fSmootherLine_grid_point(tau_gs_indx, Y, maturities_gs_indx, bw, onb)
% Performs the local linear line smoother on Y, evaluates at point m

% shortcut for the grid
g = onb.gridSpace;

% find those maturities that are closer to "m" than bw
deviations = g(tau_gs_indx) - g(maturities_gs_indx);

% if nothing close, estimate by zero
if min(abs(deviations)) >= bw
    mu_est = 0;
else
    % here perform the smoother
    
    % spatial weights for each maturity
    % Epanechnikov kernel
    W = 3/4 * ( 1-(deviations/bw).^2 );
    W( abs(deviations) >= bw ) = 0; % if deviation is greater than bw, than no weight
    
    % Gaussian kernel
%     W = normpdf(deviations/bw);
    
    % multiply by the weights corresponding to the number of available measurements
    W = W .* sum(~isnan(Y));
    
    % create the design matrix
    desM = zeros( length(maturities_gs_indx), 1 );
    desM(:,1) = 1;
%     desM(:,2) = deviations;
    
    % the average value of Y at each maturity
    Y_avg = mean(Y, 'omitnan'); 
    
    % what columns to keep (drop those where zero observations)
    keep = ~isnan(Y_avg);
    desM = desM(keep, :);
    W = diag(W(keep));
    Y_avg = Y_avg(keep);
    
    % perform the smoother
    beta=pinv(desM'*W*desM) * (desM'*W* Y_avg' );
    mu_est = beta(1);    
    
end



end