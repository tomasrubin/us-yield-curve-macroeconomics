function specCrossDensity = fEstimateSpecCrossDensity_holdout_v5(censor,onb,response,response_mean,mu_est_inSpace,bw_F,nGridFreq,holdout)

% holdout = percentage of how much data at the end of the response process is not used for the training
% typically, holdout = 0.8
holdout_number = floor(censor.nGridTime * holdout); % the number of points in the holdout partition
holdout_start = censor.nGridTime - holdout_number + 1; % the index of the first datum in the holdout partition


% pre-alocate the struct
specCrossDensity = [];
specCrossDensity.nGridFreq = nGridFreq;
specCrossDensity.gridFreq = linspace(-pi,pi,specCrossDensity.nGridFreq);

% the selection rule for Q_Bartlett is the same as for the process X only
% specCrossDensity.q_bartlett = floor(censor.nGridTime^(1/3) * (sum(censor.nGrid)/censor.nGridTime )^(1/4));
specCrossDensity.q_bartlett = floor( 2* censor.nGridTime^(1/3));
numOfLags = specCrossDensity.q_bartlett; % mind the fact that there is also lag 0; so numOfLags=3 corresponds to lags 0,1,2


%%% extract the raw covariances
% calculate the number of raw covariances - for validation
disp('crossDensity: calculate the number of raw covariances - for each lag')        
nRawCovs = zeros(1,2*numOfLags-1); 

% it's always Z_{t+lag} with X_t, but only if t+lag <(sharply) holdout_start
for lag_real = (-numOfLags+1):(numOfLags-1)
    lag_index = lag_real + numOfLags;
    lower_bound = max(1,1-lag_real);
    upper_bound = min(censor.nGridTime-lag_real,holdout_start-1); % holdout bound
    nRawCovs(lag_index) = sum( censor.nGrid( lower_bound:upper_bound) );
end

% get raw covariances (without the diagonal)
disp('crossDensity: get raw covariances (without the diagonal)')
X_jj = cell(1,2*numOfLags-1);
X_jj_onGrid = cell(1,2*numOfLags-1);
Y = cell(1,2*numOfLags-1);
for lag_real = (-numOfLags+1):(numOfLags-1)
    lag_index = lag_real + numOfLags;
    iRawCov = 1;
    X_jj{lag_index} = zeros(1, nRawCovs(lag_index));
    X_jj_onGrid{lag_index} = zeros(1, nRawCovs(lag_index));
    Y{lag_index} = zeros(1, nRawCovs(lag_index));
    
    % what interval I go over so that Z_{t+lag_real} is defined
    lower_bound = max(1,1-lag_real);
    upper_bound = min(censor.nGridTime-lag_real,holdout_start-1); % holdout bound
    for t = lower_bound:upper_bound
        for jj = 1:censor.nGrid(t)
            X_jj{lag_index}(iRawCov) = censor.grid{t}(jj);
            X_jj_onGrid{lag_index}(iRawCov) = onb.gridSpace(censor.grid{t}(jj));
            Y{lag_index}(iRawCov) = ( response(t+lag_real) - response_mean )*...
                (censor.data{t}(jj) - mu_est_inSpace(censor.grid{t}(jj)) ); % the raw covariance
            iRawCov = iRawCov + 1;            
        end
    end
end  


%%% calculate save_S save_Q

% here I'm going to save the tripple sum for Q's
% .. first index:  \nu = -specDensity.q_bartlett+1,...,+specDensity.q_bartlett-1
% .. second index: r:  1 = (r=0), 2 = (r=1)
% .. 3rd: evaluation on the onb.gridSmoother
save_Q = zeros(2*specCrossDensity.q_bartlett-1, 2, onb.nGridSmoother);

% here I'm saving the S
% .. 1st:  \nu = -specDensity.q_bartlett+1,...,+specDensity.q_bartlett-1
% .. 2nd index: pq:  0 = (r=0), 1 = (r=1), 2 = (r=2)
% .. 3rd: evaluation on the onb.gridSmoother
save_S = zeros(2*specCrossDensity.q_bartlett-1, 3, onb.nGridSmoother );

for lag_real=(-specCrossDensity.q_bartlett+1):(specCrossDensity.q_bartlett-1)
    disp(['calculating save Q and save S ... ',num2str(lag_real)])
    lag_index = lag_real + specCrossDensity.q_bartlett; % the index of the real lag in for save_S and save_Q
    
    % calculating Q^{(h)} and S^{(h)}
    for tau_jj = 1:onb.nGridSmoother
        tau_jj_real = onb.gridSmoother(tau_jj); % the real position of tau_kk on [0,1] for the smoother
        
        % prepare the spatial distance weights
        dev = X_jj_onGrid{lag_index} - tau_jj_real;
        dist = abs(dev);
        kern = 1/bw_F .* 3/4.* (1-( dist /bw_F).^2);
        kern(dist > bw_F) = 0;
        
        % save Q
        for r = 0:1
            save_Q(lag_index,r+1,tau_jj) = sum(1/censor.nGridTime .* ...
                Y{lag_index} .* ...
                kern .* ...
                ((dev)/bw_F).^r);
        end
        
        % save S
        for r = 0:2
            save_S(lag_index,r+1,tau_jj) = sum(1/censor.nGridTime .* ...
                kern .* ...
                ((dev)/bw_F).^r);
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%         PERFOM THE MANIPULATION TO EVALUATE THE DENSITY ON THE FREQ GRID             %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare S
S_0 = zeros(onb.nGridSmoother,1);
S_1 = zeros(onb.nGridSmoother,1);
S_2 = zeros(onb.nGridSmoother,1);

for lag=(-specCrossDensity.q_bartlett+1):(specCrossDensity.q_bartlett-1)
    lag_index = lag + specCrossDensity.q_bartlett;
    S_0 = S_0 + squeeze( save_S(lag_index,1,:) ) * (1-abs(lag)/specCrossDensity.q_bartlett);
    S_1 = S_1 + squeeze( save_S(lag_index,2,:) ) * (1-abs(lag)/specCrossDensity.q_bartlett);
    S_2 = S_2 + squeeze( save_S(lag_index,3,:) ) * (1-abs(lag)/specCrossDensity.q_bartlett);    
end


% estimate the spec density using the Bartlett's weights
specCrossDensity.smoother = zeros(specCrossDensity.nGridFreq, onb.nGridSmoother); % this is in the grid for which I smoothed the lag covariances
specCrossDensity.ONB = zeros(specCrossDensity.nGridFreq, onb.nBasis);
for k = 1:specCrossDensity.nGridFreq
    omega = specCrossDensity.gridFreq(k); % real omega: -pi to pi       
    Q_0 = zeros(onb.nGridSmoother,1);
    Q_1 = zeros(onb.nGridSmoother,1);
    for lag=(-specCrossDensity.q_bartlett+1):(specCrossDensity.q_bartlett-1)
        lag_index = lag + specCrossDensity.q_bartlett;
        Q_0 = Q_0 + squeeze(save_Q(lag_index, 1, :)) * exp(-1i*lag*omega) * (1-abs(lag)/specCrossDensity.q_bartlett);
        Q_1 = Q_1 + squeeze(save_Q(lag_index, 2, :)) * exp(-1i*lag*omega) * (1-abs(lag)/specCrossDensity.q_bartlett);
    end

    % the spectral CROSS-density
    m = 1/(2*pi) * specCrossDensity.q_bartlett * (Q_0.*S_2 - S_1.*Q_1) ./ (S_0.*S_2-S_1.^2); % mind the fact that the sum of Batlett's weights is specDensity.q_bartlett    
    specCrossDensity.smoother(k,:) = m;
    specCrossDensity.ONB(k,:) = onb.smoothingMatrix_gridSmoother2ONB*m;
end


end