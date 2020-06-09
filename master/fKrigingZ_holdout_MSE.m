function mse = fKrigingZ_holdout_MSE( response_est, response , holdout)

% holdout = percentage of how much data at the end of the response process is not used for the training
% typically, holdout = 0.8
nGridTime = length(response);
holdout_number = floor(nGridTime * holdout); % the number of points in the holdout partition
holdout_start = nGridTime - holdout_number + 1; % the index of the first datum in the holdout partition

mse = mean( (response(holdout_start:end)-response_est(holdout_start:end)).^2 );

end