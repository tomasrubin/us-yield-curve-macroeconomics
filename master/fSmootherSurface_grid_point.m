function R_est = fSmootherSurface_grid_point(tau1_gs_indx, tau2_gs_indx, raws_G, raws_n, maturities_gs_indx, bw_1, bw_2, onb)
% Performs the local linear line smoother on Y, evaluates at point m

% shortcut for the grid
g = onb.gridSpace;

% find those maturities that are closer to "m" than bw
deviations_1 = g(tau1_gs_indx) - g(maturities_gs_indx);
deviations_2 = g(tau2_gs_indx) - g(maturities_gs_indx);

deviation_matrix_1 = repmat(deviations_1,  length(maturities_gs_indx), 1);
deviation_matrix_2 = repmat(deviations_2', 1, length(maturities_gs_indx));

% if nothing close, estimate by zero
if (min(abs(deviations_1)) >= bw_1) || (min(abs(deviations_2)) >= bw_2)
    R_est = 0;
else
    % here perform the smoother
    
    % spatial weights for each maturity
    % Epanechnikov kernel
    W = 3/4 * ( 1-(deviation_matrix_1(:)/bw_1).^2 ) .* ( 1-(deviation_matrix_2(:)/bw_2).^2 );
    W( abs(deviation_matrix_1(:)) >= bw_1 ) = 0; % if deviation is greater than bw, than no weight
    W( abs(deviation_matrix_2(:)) >= bw_2 ) = 0; % if deviation is greater than bw, than no weight
    
    % Gaussian kernel
%     W = normpdf(deviation_matrix_1(:)/bw_1) .* normpdf(deviation_matrix_2(:)/bw_2);
    
    % multiply by the weights corresponding to the number of available measurements
    W = W .* sum(raws_n(:));
    
    % create the design matrix
    desM = zeros( length(maturities_gs_indx).^2, 2 );
    desM(:,1) = 1;
    desM(:,2) = deviation_matrix_1(:);
    desM(:,3) = deviation_matrix_2(:);
        
    % what columns to keep (drop those where zero observations)
    keep = ( raws_n(:) ~= 0 );
    desM = desM(keep, :);
    W = diag(W(keep));
    raws_G = raws_G(keep);
    
    % perform the smoother
    beta=pinv(desM'*W*desM) * (desM'*W* raws_G );
    R_est = beta(1);    
    
end



end