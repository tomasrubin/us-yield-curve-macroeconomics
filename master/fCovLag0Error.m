function covLag0_error = fCovLag0Error( onb, covsLags_true, covLag0_ONB )

covLag0_error = [];

% it ONB
covLag0_true = covsLags_true.lag0;

%% cov lag 0 error
m = onb.onbMatrix * (covLag0_true - covLag0_ONB) * onb.onbMatrix';                
m2 = onb.onbMatrix*covLag0_true*onb.onbMatrix';
covLag0_error.mse = norm(m, 'fro')^2 / onb.nGridSpace^2;
covLag0_error.rmse = norm(m, 'fro')^2 / norm(m2, 'fro')^2;

end