function [bw_mu_opt,mu_est_smoother_opt] = fKCVforMU(censor,onb)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% pre-calculate some quantities

% data counts for training
nData_train = nan(censor.cv_K,1);
for kk = 1:censor.cv_K, nData_train(kk) = sum(censor.cv_batch_counts([1:kk-1 kk+1:end])); end

% data counts for validation = identical to what is in each bach
nData_valid = censor.cv_batch_counts;

% prepare the data pool for K-folds
Xs_train = cell(censor.cv_K,1); % here will be the values on [0,1] - for training - all except partition k
Xs_valid = cell(censor.cv_K,1); % for validation - partition k only
Xg_train = cell(censor.cv_K,1); % intex of the grid
Xg_valid = cell(censor.cv_K,1);
Y_train = cell(censor.cv_K,1); % value of the function
Y_valid = cell(censor.cv_K,1);

% save all data - it will be used for the final fit on all data
nData_all = sum(censor.cv_batch_counts);
Xs_all = nan( nData_all,1 );
Xg_all = nan( nData_all,1 );
Y_all = nan( nData_all,1 );


% cycle for each partition
for kk = 1:censor.cv_K    
    % prepare data structures
    Xs_train{kk} = nan( nData_train(kk) , 1);
    Xg_train{kk} = nan( nData_train(kk) , 1);
    Y_train{kk} = nan( nData_train(kk) , 1);
    Xs_valid{kk} = nan( nData_valid(kk) , 1);
    Xg_valid{kk} = nan( nData_valid(kk) , 1);
    Y_valid{kk} = nan( nData_valid(kk) , 1);
    
    % counters of current data points
    iData_train = 1;
    iData_valid = 1;
    iData_all = 1;
    
    % fill-in the data
    for tt = 1:censor.nGridTime
        for iii = 1:censor.nGrid(tt)
            if censor.cv_batch{tt}(iii) ~= kk
                % TRAINING SET : include iff it's not "kk"
                Xs_train{kk}(iData_train) = onb.gridSpace(censor.grid{tt}(iii));
                Xg_train{kk}(iData_train) = censor.grid{tt}(iii);
                Y_train{kk}(iData_train) = censor.data{tt}(iii);
                iData_train = iData_train + 1;
            else
                % VALIDATION SET : include iff it's "kk"
                Xs_valid{kk}(iData_valid) = onb.gridSpace(censor.grid{tt}(iii));
                Xg_valid{kk}(iData_valid) = censor.grid{tt}(iii);
                Y_valid{kk}(iData_valid) = censor.data{tt}(iii);
                iData_valid = iData_valid + 1;
            end
            % ALL CASES
            Xs_all(iData_all) = onb.gridSpace(censor.grid{tt}(iii));
            Xg_all(iData_all) = censor.grid{tt}(iii);
            Y_all(iData_all) = censor.data{tt}(iii);
            iData_all = iData_all + 1;
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Baysian optimization

% perform the optimization
f = @(x)sum_of_squares(x.bw_mu); % I don't know why but I have to create the anonymous function handle
bw_mu = optimizableVariable('bw_mu',[0.02,.3]);
results = bayesopt(f, bw_mu, 'PlotFcn', [], 'MaxObjectiveEvaluations', 10, 'IsObjectiveDeterministic', true);

% return the optimal BW and its estimated mu
bw_mu_opt = results.XAtMinObjective.bw_mu;


%% fit again on full data
mu_est_smoother_opt = fSmootherLine( Xs_all, Y_all, bw_mu_opt, 1,onb.gridSmoother  );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to calculate the CV score for a given candidate bandwidth
    function mse = sum_of_squares(bw_mu)

        mse = 0; % the mean square error
        for k = 1:censor.cv_K
            
            % perform the line smoother - on the training data partition except k
            mu_est = fSmootherLine(Xs_train{k},Y_train{k},bw_mu,1,onb.gridSmoother);
            mu_est_inSpace = onb.onbMatrix*onb.smoothingMatrix_gridSmoother2ONB* mu_est;
            
            % calculate the sum of squares - on validation partition k
            mse = mse + 1/censor.cv_K * mean(( mu_est_inSpace(Xg_valid{k}) - Y_valid{k} ).^2);
        end
        
    end

end