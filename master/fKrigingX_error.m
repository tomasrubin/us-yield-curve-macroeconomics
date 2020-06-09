function error = fKrigingX_error( onb, dataFull, krigingX, covLag0, checkCoverage )
% This function checks the kriging error
% onb
% dataFull ... the true latent functional data (inSpace)
% krigingX ... struct created while kriging (in ONB)
% checkCoverage ... 1 = also check the coverage of pointwise and simul
%                       confidence bands
% covLag0 ... true cov lag 0 for RMSE


% pre-alocate space
error = [];
error.mse = 0;
error.rmse = 0;
error.coverage_point = 0;
error.coverage_simul = 0;

[nGridTime,~] = size(dataFull);
for t = 1:nGridTime
    
    % kriging error
    error.mse  = error.mse  + 1/nGridTime* mean( (dataFull(t,:)' - onb.onbMatrix*real(krigingX.est(t,:))').^2);
    error.rmse = error.rmse + 1/nGridTime* mean( (dataFull(t,:)' - onb.onbMatrix*real(krigingX.est(t,:))').^2) ./ trace(covLag0);
    
    % confidence bands coverage check
    if checkCoverage
        
        % just short-cuts
        post_mu = real(krigingX.est(t,:));
        post_cov = real(squeeze(krigingX.var(t,:,:)));
        
        % pointwise band
        lower_point = onb.onbMatrix*post_mu' - 1.96*sqrt(abs(diag( onb.onbMatrix*post_cov*onb.onbMatrix' )));
        upper_point = onb.onbMatrix*post_mu' + 1.96*sqrt(abs(diag( onb.onbMatrix*post_cov*onb.onbMatrix' )));    
        
        % simultaneaous band
        z = fSCBDegras( post_cov, onb.onbMatrix );
        lower_simul = onb.onbMatrix*post_mu' - z*sqrt(abs(diag( onb.onbMatrix*post_cov*onb.onbMatrix')));        
        upper_simul = onb.onbMatrix*post_mu' + z*sqrt(abs(diag( onb.onbMatrix*post_cov*onb.onbMatrix')));
        
        % quantify the coverage
        error.coverage_point = error.coverage_point + 1/nGridTime * mean( (dataFull(t,:)' > lower_point) & (dataFull(t,:)' < upper_point) );
        error.coverage_simul = error.coverage_point + 1/nGridTime * prod( (dataFull(t,:)' > lower_simul) & (dataFull(t,:)' < upper_simul) );
        
    end
    
end






end