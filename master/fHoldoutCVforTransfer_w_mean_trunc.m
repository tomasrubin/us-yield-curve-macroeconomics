function [threshold_value_opt] = fHoldoutCVforTransfer_w_mean_trunc(threshold_value_interval, censor, onb, CV_on_proportion,specDensity,specCrossDensity,krigingX_dyna_padded,response,dataFull_padding,mu_ONB,response_mean)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baysian optimization

% perform the optimization
f = @(x)holdout_mse(x.log_coefficient); % I don't know why but I have to create the anonymous function handle
log_coefficient = optimizableVariable('log_coefficient', log(threshold_value_interval) );
results = bayesopt(f, log_coefficient, 'PlotFcn', [], 'MaxObjectiveEvaluations', 20, 'IsObjectiveDeterministic', true);
% results = bayesopt(f, log_coefficient, 'MaxObjectiveEvaluations', 20, 'IsObjectiveDeterministic', true);


% return the optimal BW and its estimated mu
log_coefficient_opt = results.XAtMinObjective.log_coefficient;
threshold_value_opt = exp(log_coefficient_opt) * censor.nGridTime.^(-1/3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to calculate the CV score for a given candidate coefficient

    function x_mse = holdout_mse(x_log_coefficient)        

        % fit the transfer functions
        x_threshold_value = exp(x_log_coefficient) * censor.nGridTime.^(-1/3);
        x_transfer = fEstimateTransfer_threshold_value(onb,specDensity,specCrossDensity, x_threshold_value);

        % estimation on the TRAINING sample
        x_response_est = fKrigingZ_nonzeromean(krigingX_dyna_padded, x_transfer, mu_ONB, response_mean, dataFull_padding);
        
        % MSE on the HOLDOUT
        x_mse = fKrigingZ_holdout_MSE_NaNready( x_response_est, response', CV_on_proportion);
        
    end

end