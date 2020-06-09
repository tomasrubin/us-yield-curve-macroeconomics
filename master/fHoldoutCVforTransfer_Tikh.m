function [tikhonovParameter_opt] = fHoldoutCVforTransfer_Tikh(tikhonovParameter_interval, censor, onb, holdout,specDensity,specCrossDensity,krigingX_dyna_padded,response,dataFull_padding)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baysian optimization

% perform the optimization
f = @(x)holdout_mse(x.log_coefficient); % I don't know why but I have to create the anonymous function handle
log_coefficient = optimizableVariable('log_coefficient',log(tikhonovParameter_interval) );
results = bayesopt(f, log_coefficient, 'PlotFcn', [], 'MaxObjectiveEvaluations', 20, 'IsObjectiveDeterministic', true);
% results = bayesopt(f, log_coefficient, 'MaxObjectiveEvaluations', 40, 'IsObjectiveDeterministic', true);

% return the optimal BW and its estimated mu
log_coefficient_opt = results.XAtMinObjective.log_coefficient;
tikhonovParameter_opt = exp(log_coefficient_opt) * censor.nGridTime ^(-1/3) * (sum(censor.nGrid)/censor.nGridTime)^(-1/3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to calculate the CV score for a given candidate coefficient

    function x_mse = holdout_mse(x_log_coefficient)        

        % fit the transfer functions
        x_tikhonovParameter = exp(x_log_coefficient) * censor.nGridTime ^(-1/3) * (sum(censor.nGrid)/censor.nGridTime)^(-1/3);        
        x_transfer = fEstimateTransfer_Tikhonov(onb,specDensity,specCrossDensity, x_tikhonovParameter);

        % estimation on the TRAINING sample
        x_response_est = fKrigingZ(krigingX_dyna_padded,x_transfer,dataFull_padding);
        
        % MSE on the HOLDOUT
        x_mse = fKrigingZ_holdout_MSE( x_response_est, response', holdout);
        
    end

end