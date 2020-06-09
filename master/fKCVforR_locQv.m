function [bw_r_opt, bw_v_opt, covLag0_smoother, covLag0_ridge_smoother, sigma2_est_loc_qv, sigma2_est_loc_lin] = fKCVforR_locQv(censor,onb,mu_est_inSpace)


%% precalculate some stuff

% here go all raw covariances excluding k - and off diagonal - for smoothing the kernel
nRawCovs_all_but_k = zeros(censor.cv_K,1); % the number
X1_all_but_k = cell(censor.cv_K,1); % the first covariate
X2_all_but_k = cell(censor.cv_K,1); % the second covariate
Y_all_but_k = cell(censor.cv_K,1); % the value of the raw covariates product

% here go all raw covariances excluding "k" on the diagonal - for smoothing
nRawCovs_diag_but_k = zeros(censor.cv_K,1); % the number
XX_diag_but_k = cell(censor.cv_K,1); % the diagonal coordinate
Y_diag_but_k = cell(censor.cv_K,1); % the value of the raw covariates product

% here go raw covariances from "k" - for CV check (both diag and off diag)
nRawCovs_diag_k_4_cv = zeros(censor.cv_K,1);
XX_ind_diag_k_4_cv = cell(censor.cv_K,1); % the first covariate
Y_diag_k_4_cv = cell(censor.cv_K,1); % the value of the raw covariates product
nRawCovs_off_k_4_cv = zeros(censor.cv_K,1);
X1_ind_off_k_4_cv = cell(censor.cv_K,1); % the first covariate
X2_ind_off_k_4_cv = cell(censor.cv_K,1); % the second covariate
Y_off_k_4_cv = cell(censor.cv_K,1); % the value of the raw covariates product

% cycle over all partitions
for cv_kk = 1:censor.cv_K

    % calculate the number of raw covariances    
    for tt = 1:censor.nGridTime
        for iii = 1:censor.nGrid(tt)
            for jjj = 1:censor.nGrid(tt)
                
                % off-diagonal and not in the "cv_kk"
                if (iii~=jjj) && (censor.cv_batch{tt}(iii) ~= cv_kk) && (censor.cv_batch{tt}(jjj) ~= cv_kk)
                    nRawCovs_all_but_k(cv_kk) = nRawCovs_all_but_k(cv_kk)+1;
                end
                % on-diagonal and not in the "cv_kk"
                if (iii==jjj) && (censor.cv_batch{tt}(iii) ~= cv_kk) && (censor.cv_batch{tt}(jjj) ~= cv_kk)
                    nRawCovs_diag_but_k(cv_kk) = nRawCovs_diag_but_k(cv_kk)+1;
                end
                % in cv_kk --- for validation, diag
                if (iii==jjj) && (censor.cv_batch{tt}(iii) == cv_kk) && (censor.cv_batch{tt}(jjj) == cv_kk)
                    nRawCovs_diag_k_4_cv(cv_kk) = nRawCovs_diag_k_4_cv(cv_kk)+1;
                end
                % in cv_kk --- for validation, off diag
                if (iii~=jjj) && (censor.cv_batch{tt}(iii) == cv_kk) && (censor.cv_batch{tt}(jjj) == cv_kk)
                    nRawCovs_off_k_4_cv(cv_kk) = nRawCovs_off_k_4_cv(cv_kk)+1;
                end
            end
        end
    end

    % all raw covariances but kk - off diagonal
    X1_all_but_k{cv_kk} = zeros(nRawCovs_all_but_k(cv_kk),1); % the "covariates"
    X2_all_but_k{cv_kk} = zeros(nRawCovs_all_but_k(cv_kk),1); % the "covariates"
    Y_all_but_k{cv_kk} = zeros(nRawCovs_all_but_k(cv_kk),1); % the value of the raw covariates product
    
    % all raw covariances but kk - on diagonal
    XX_diag_but_k{cv_kk} = zeros(nRawCovs_diag_but_k(cv_kk),1); % the "covariates"
    Y_diag_but_k{cv_kk} = zeros(nRawCovs_diag_but_k(cv_kk),1); % the value of the raw covariates product
    
    % in kk
    XX_ind_diag_k_4_cv{cv_kk} = zeros(nRawCovs_diag_k_4_cv(cv_kk),1); % the "covariates"
    Y_diag_k_4_cv{cv_kk} = zeros(nRawCovs_diag_k_4_cv(cv_kk),1); % the value of the raw covariates product
    
    X1_ind_off_k_4_cv{cv_kk} = zeros(nRawCovs_off_k_4_cv(cv_kk),1); % the "covariates"
    X2_ind_off_k_4_cv{cv_kk} = zeros(nRawCovs_off_k_4_cv(cv_kk),1); % the "covariates"
    Y_off_k_4_cv{cv_kk} = zeros(nRawCovs_off_k_4_cv(cv_kk),1); % the value of the raw covariates product
    
    
    iRawCov_all_but_k = 1; % the index of the raw cow now being filled
    iRawCov_diag_but_k = 1;
    iRawCov_diag_k_4_cv = 1;
    iRawCov_off_k_4_cv = 1;
    for tt = 1:censor.nGridTime
        for iii = 1:censor.nGrid(tt)
            for jjj = 1:censor.nGrid(tt)
                % off-diagonal and not in the "cv_kk"
                if (iii~=jjj) && (censor.cv_batch{tt}(iii) ~= cv_kk) && (censor.cv_batch{tt}(jjj) ~= cv_kk)
                    X1_all_but_k{cv_kk}(iRawCov_all_but_k,1) = onb.gridSpace(censor.grid{tt}(iii));
                    X2_all_but_k{cv_kk}(iRawCov_all_but_k,1) = onb.gridSpace(censor.grid{tt}(jjj));
                    Y_all_but_k{cv_kk}(iRawCov_all_but_k) = (censor.data{tt}(iii) - mu_est_inSpace(censor.grid{tt}(iii)))*...
                        (censor.data{tt}(jjj) - mu_est_inSpace(censor.grid{tt}(jjj)));
                    iRawCov_all_but_k = iRawCov_all_but_k+1;
                end
                % on-diagonal and not in the "cv_kk"
                if (iii==jjj) && (censor.cv_batch{tt}(iii) ~= cv_kk) && (censor.cv_batch{tt}(jjj) ~= cv_kk)
                    XX_diag_but_k{cv_kk}(iRawCov_diag_but_k,1) = onb.gridSpace(censor.grid{tt}(iii));
                    Y_diag_but_k{cv_kk}(iRawCov_diag_but_k) = (censor.data{tt}(iii) - mu_est_inSpace(censor.grid{tt}(iii)))*...
                        (censor.data{tt}(jjj) - mu_est_inSpace(censor.grid{tt}(jjj)));
                    iRawCov_diag_but_k = iRawCov_diag_but_k+1;
                end                
                % in cv_kk --- for validation, diag
                if (iii==jjj) && (censor.cv_batch{tt}(iii) == cv_kk) && (censor.cv_batch{tt}(jjj) == cv_kk)
                    XX_ind_diag_k_4_cv{cv_kk}(iRawCov_diag_k_4_cv,1) = censor.grid{tt}(iii);
                    Y_diag_k_4_cv{cv_kk}(iRawCov_diag_k_4_cv) = (censor.data{tt}(iii) - mu_est_inSpace(censor.grid{tt}(iii)))*...
                        (censor.data{tt}(jjj) - mu_est_inSpace(censor.grid{tt}(jjj)));
                    iRawCov_diag_k_4_cv = iRawCov_diag_k_4_cv+1;
                end
                % in cv_kk --- for validation, off diag
                if (iii~=jjj) && (censor.cv_batch{tt}(iii) == cv_kk) && (censor.cv_batch{tt}(jjj) == cv_kk)
                    X1_ind_off_k_4_cv{cv_kk}(iRawCov_off_k_4_cv,1) = censor.grid{tt}(iii);
                    X2_ind_off_k_4_cv{cv_kk}(iRawCov_off_k_4_cv,1) = censor.grid{tt}(jjj);
                    Y_off_k_4_cv{cv_kk}(iRawCov_off_k_4_cv) = (censor.data{tt}(iii) - mu_est_inSpace(censor.grid{tt}(iii)))*...
                        (censor.data{tt}(jjj) - mu_est_inSpace(censor.grid{tt}(jjj)));
                    iRawCov_off_k_4_cv = iRawCov_off_k_4_cv+1;
                end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here go all raw covariances

% off diag
nRawCovs_all = sum(censor.nGrid.^2 - censor.nGrid);
X1_all = zeros(nRawCovs_all,1);
X2_all = zeros(nRawCovs_all,1);
Y_all = zeros(nRawCovs_all,1);
% on diag
nRawCovs_diag_all = sum(censor.nGrid);
XX_diag_all = zeros(nRawCovs_diag_all,1);
Y_diag_all = zeros(nRawCovs_diag_all,1);

% running indexes
iRawCov_all = 1;
iRawCov_diag = 1;

for tt = 1:censor.nGridTime
    for iii = 1:censor.nGrid(tt)
        for jjj = 1:censor.nGrid(tt)
            
            if iii~=jjj % off-diagonal
                X1_all(iRawCov_all) = onb.gridSpace(censor.grid{tt}(iii));
                X2_all(iRawCov_all) = onb.gridSpace(censor.grid{tt}(jjj)); 
                Y_all(iRawCov_all) = (censor.data{tt}(iii) - mu_est_inSpace(censor.grid{tt}(iii)))*...
                        (censor.data{tt}(jjj) - mu_est_inSpace(censor.grid{tt}(jjj)));
                iRawCov_all = iRawCov_all + 1;
            else % on-diagonal
                XX_diag_all(iRawCov_diag) = onb.gridSpace(censor.grid{tt}(iii));
                Y_diag_all(iRawCov_diag) = (censor.data{tt}(iii) - mu_est_inSpace(censor.grid{tt}(iii)))*...
                        (censor.data{tt}(jjj) - mu_est_inSpace(censor.grid{tt}(jjj)));
                iRawCov_diag = iRawCov_diag + 1;
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BayesOpt - sum_of_squares_off_diag
f1 = @(x)sum_of_squares_off_diag(x.bw_r);
bw_r = optimizableVariable('bw_r',[0.01,0.5]);

% return the optimal BW
results1 = bayesopt(f1, bw_r, 'PlotFcn', [],'verbose',1, 'MaxObjectiveEvaluations', 20, 'IsObjectiveDeterministic', true);
if isnan(results1.MinObjective)% too sparse design
    bw_r_opt = 1;
else
    bw_r_opt = results1.XAtMinObjective.bw_r;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BayesOpt - sum_of_squares_on_diag
f2 = @(x)sum_of_squares_on_diag(x.bw_v);
bw_v = optimizableVariable('bw_v',[0.005,0.9]);

% return the optimal BW
results2 = bayesopt(f2, bw_v, 'PlotFcn', [],'verbose',1, 'MaxObjectiveEvaluations', 20, 'IsObjectiveDeterministic', true);
if isnan(results2.MinObjective) % too sparse design
    bw_v_opt = 1;
else
    bw_v_opt = results2.XAtMinObjective.bw_v;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% re-estimate quantities from all data

% cov lag 0
covLag0_smoother = fSmootherSurface( [X1_all,X2_all], Y_all, bw_r_opt, onb.gridSmoother, 1 );

% smothDiag
covLag0_ridge_smoother = fSmootherLine( XX_diag_all, Y_diag_all, bw_v_opt, 1, onb.gridSmoother );

%% estimate the sigma
% local quadratic estimate around diagonal

% precalculate the loc-quadratic fit on the diagonal
covLag0_diag_qv = zeros(onb.nGridSmoother,1);
for iii = 1:onb.nGridSmoother
    % now we're filling up the point covLag0(iii,iii), i.e. [onb.gridSpace(iii), onb.gridSpace(iii)]                
    % choose only those indeces where I'm less then "bw" far away
    indxo = find( abs(X1_all - onb.gridSmoother(iii)) <= bw_r_opt + 10^(-6) & abs( X2_all - onb.gridSmoother(iii)) <= bw_r_opt + 10^(-6) );
    dist1o = X1_all(indxo) - onb.gridSmoother(iii);        
    dist2o = X2_all(indxo) - onb.gridSmoother(iii); 
    diag_closesto = (X1_all(indxo)+X2_all(indxo))/2;
    dist_to_diag = sqrt( (X1_all(indxo) - diag_closesto).^2 + (X2_all(indxo) - diag_closesto).^2 );    
    dist_to_diag( X2_all(indxo) < X1_all(indxo) ) = -dist_to_diag( X2_all(indxo) < X1_all(indxo) ); % make it negative when I'm in the lower-right triangle    
    % define the weights - it's the same
    Wo = spdiags( 9/16*(1-(dist1o/bw_r_opt).^2).*(1-(dist2o/bw_r_opt).^2), 0, length(indxo),length(indxo) );
%     W = diag( 9/16*(1-(dist1/bw_r).^2).*(1-(dist2/bw_r).^2));
    % get the design matrix - only the distance from the diagonal matters
    desMo = zeros( length(indxo), 3 );
    desMo(:,1) = 1;
    desMo(:,2) = dist_to_diag;
    desMo(:,3) = dist_to_diag.^2;
    betao=pinv(desMo'*Wo*desMo) * (desMo'*Wo* Y_all(indxo) );
    covLag0_diag_qv(iii) = betao(1);
end

% express in Space
covLag0_est_inSpace = onb.onbMatrix * onb.smoothingMatrix_gridSmoother2ONB * covLag0_smoother *onb.smoothingMatrix_gridSmoother2ONB'*onb.onbMatrix';
covLag0_ridge_inSpace = onb.onbMatrix * onb.smoothingMatrix_gridSmoother2ONB * covLag0_ridge_smoother;
covLag0_diag_qv_inSpace = onb.onbMatrix * onb.smoothingMatrix_gridSmoother2ONB * covLag0_diag_qv;

% local quadratic smoothing on the diagonal
sigma2_est_loc_qv = mean(covLag0_ridge_inSpace-covLag0_diag_qv_inSpace);

% local linear on the diagonal - uses the estimate from the surface smoother for covLag0
sigma2_est_loc_lin = mean( covLag0_ridge_inSpace - diag(covLag0_est_inSpace) );











%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to calculate the CV score for a given candidate bandwidth - OFF DIAGONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function mse = sum_of_squares_off_diag(bw_r)
           
        mse = 0; % the mean square error
        for k = 1:censor.cv_K
            
            % surface smoother
            covLag0 = fSmootherSurface( [X1_all_but_k{k},X2_all_but_k{k}], Y_all_but_k{k},bw_r,onb.gridSmoother,1 );
            covLag0_onb = onb.smoothingMatrix_gridSmoother2ONB*covLag0*onb.smoothingMatrix_gridSmoother2ONB'; % convert into the ONB
            covLag0_inSpace = onb.onbMatrix * covLag0_onb * onb.onbMatrix'; % without the diagonal augmentation
            
            % prepare the fitted values
            Y_estim_off = nan( size(Y_off_k_4_cv{k}) );
            for  jj = 1:length(Y_estim_off)
                Y_estim_off(jj) = covLag0_inSpace( X1_ind_off_k_4_cv{k}(jj), X2_ind_off_k_4_cv{k}(jj));
            end
            
            % calculate the CV score
            mse = mse + 1/censor.cv_K * ( mean((Y_estim_off - Y_off_k_4_cv{k}).^2) );
            
        end % of forcycle for cross-validation
        
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to calculate the CV score for a given candidate bandwidth - ON DIAGONAL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mse = sum_of_squares_on_diag(bw_v)
                
        mse = 0; % the mean square error
        for k = 1:censor.cv_K
            
            % line smoother
            smooth_diag = fSmootherLine( XX_diag_but_k{k}, Y_diag_but_k{k}, bw_v, 1, onb.gridSmoother );
            smooth_diag_onb = onb.smoothingMatrix_gridSmoother2ONB * smooth_diag;
            smooth_diag_inSpace = onb.onbMatrix*smooth_diag_onb;
            
            % calculate the CV score
            Y_estim_diag = smooth_diag_inSpace( XX_ind_diag_k_4_cv{k});            
            mse = mse + 1/censor.cv_K * ( mean((Y_estim_diag - Y_diag_k_4_cv{k}).^2) );
            
        end % of forcycle for cross-validation        
    end

end