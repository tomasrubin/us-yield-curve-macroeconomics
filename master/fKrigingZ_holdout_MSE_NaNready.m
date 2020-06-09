function mse = fKrigingZ_holdout_MSE_NaNready( response_est, response , error_calculate_on_proportion)

% series length
nGridTime = length(response);

% where to start
CV_start = ceil( nGridTime * error_calculate_on_proportion(1) + 0.001);

% where to end
CV_end = ceil( nGridTime * error_calculate_on_proportion(2));

% disp("fKrigingZ_holdout_MSE_NaNready: MSE on "+num2str(CV_start)+"..."+num2str(CV_end))

x = response(CV_start:CV_end)-response_est(CV_start:CV_end);
mse = mean( (x).^2, 'omitnan' );

end