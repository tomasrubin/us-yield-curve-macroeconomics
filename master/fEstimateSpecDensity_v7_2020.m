function specDensity = fEstimateSpecDensity_v7_2020(censor,onb,mu_est_inSpace,bw_F,nGridFreq)
% this version has different rule for Bartletts span parameter

specDensity = [];
specDensity.nGridFreq = nGridFreq;
specDensity.gridFreq = linspace(-pi,pi,specDensity.nGridFreq);

% selection rule for specDensity.q_bartlett
specDensity.q_bartlett = floor(censor.nGridTime^(1/3) * (sum(censor.nGrid)/censor.nGridTime )^(1/4));
% specDensity.q_bartlett = floor( 2* censor.nGridTime^(1/3));
numOfLags = specDensity.q_bartlett; % mind the fact that there is also lag 0; so numOfLags=3 corresponds to lags 0,1,2

%%% extract the raw covariances
% calculate the number of raw covariances - for validation
disp('calculate the number of raw covariances - for each lag')        
nRawCovs = zeros(1,numOfLags); 

for t = 1:censor.nGridTime
    for s = t:min(censor.nGridTime, t+numOfLags-1)
        lag = s-t;
        for ii = 1:censor.nGrid(t)
            for jj = 1:censor.nGrid(s)
                if (ii ~= jj) || (t ~= s) % calculate only off-diagonal if lag-0
                    nRawCovs(lag+1) = nRawCovs(lag+1)+1;
                end
            end
        end
    end
end        

% get raw covariances (without the diagonal)
disp('get raw covariances (without the diagonal)')
X_ii = cell(1,numOfLags);
X_jj = cell(1,numOfLags);
X_ii_onGrid = cell(1,numOfLags);        
X_jj_onGrid = cell(1,numOfLags);
Y = cell(1,numOfLags);
for lag = 0:(numOfLags-1)
    iRawCov = 1;
    X_ii{lag+1} = zeros(1, nRawCovs(lag+1));
    X_jj{lag+1} = zeros(1, nRawCovs(lag+1));
    X_ii_onGrid{lag+1} = zeros(1, nRawCovs(lag+1));
    X_jj_onGrid{lag+1} = zeros(1, nRawCovs(lag+1));
    Y{lag+1} = zeros(1, nRawCovs(lag+1));
    for t = 1:(censor.nGridTime-lag)
        for ii = 1:censor.nGrid(t+lag)
            for jj = 1:censor.nGrid(t)
                if (ii ~= jj) || (lag ~= 0) % calculate only off-diagonal if lag-0
                    X_ii{lag+1}(iRawCov) = censor.grid{t+lag}(ii);                            
                    X_jj{lag+1}(iRawCov) = censor.grid{t}(jj);
                    X_ii_onGrid{lag+1}(iRawCov) = onb.gridSpace(censor.grid{t+lag}(ii));
                    X_jj_onGrid{lag+1}(iRawCov) = onb.gridSpace(censor.grid{t}(jj));
                    Y{lag+1}(iRawCov) = (censor.data{t+lag}(ii)-mu_est_inSpace(censor.grid{t+lag}(ii)) )*...
                        (censor.data{t}(jj) - mu_est_inSpace(censor.grid{t}(jj)) ); % the raw covariance
                    iRawCov = iRawCov + 1;
                end
            end
        end                
    end
end    


%%% estimate the average number of measurements per curve
bar_n = mean( censor.nGrid ); % average number of points per curve
bar_n2 = mean( censor.nGrid .^2 ); % average number squared per curve


%%% calculate save_S save_Q

% here I'm going to save the tripple sum for Q's
% .. first index:  \nu = -specDensity.q_bartlett+1,...,+specDensity.q_bartlett-1
% .. second index: pq:  1 = (p=0,q=0), 2 = (p=1,q=0), 3 = (p=0,q=1)
% .. 3rd & 4th : evaluation on the onb.gridSmoother
save_Q = zeros(2*specDensity.q_bartlett-1, 3, onb.nGridSmoother, onb.nGridSmoother);

% here I'm saving the S
% .. 1st:  \nu = -specDensity.q_bartlett+1,...,+specDensity.q_bartlett-1
% .. 2nd index: pq:  1 = (p=0,q=0), 2 = (p=1,q=0), 3 = (p=0,q=1), 4 = (p=1,q=1), 5 = (p=2,q=0), 6 = (p=0,q=2)
% .. 3rd & 4th : evaluation on the onb.gridSmoother
save_S = zeros(2*specDensity.q_bartlett-1, 6, onb.nGridSmoother, onb.nGridSmoother);

for lag=0:(specDensity.q_bartlett-1) 
    disp(['calculating save Q and save S ... ',num2str(lag)])            
    nu_i = lag + specDensity.q_bartlett; % the index of the real lag in for save_S and save_Q

    % calculating Q^{(h)} and S^{(h)}
    for tau_ii = 1:onb.nGridSmoother
        tau_ii_real = onb.gridSmoother(tau_ii); % the real position of tau_jj on [0,1] for the smoother
        for tau_jj = 1:onb.nGridSmoother                            
            tau_jj_real = onb.gridSmoother(tau_jj); % the real position of tau_kk on [0,1] for the smoother

            dev1 = X_ii_onGrid{lag+1} - tau_ii_real;
            dev2 = X_jj_onGrid{lag+1} - tau_jj_real;
            dist1 = abs(dev1);
            dist2 = abs(dev2);
            kern1 = 1/bw_F .* 3/4.* (1-( dist1 /bw_F).^2);                    
            kern2 = 1/bw_F .* 3/4.* (1-( dist2 /bw_F).^2);
            kern1(dist1 > bw_F) = 0;
            kern2(dist2 > bw_F) = 0;
            % save Q
            for pq = 1:3
                switch pq
                    case 1, p=0; q=0;
                    case 2, p=1; q=0;                                   
                    case 3, p=0; q=1;                                
                end
                ss = sum(...
                    Y{lag+1} .* ...
                    kern1 .* ...
                    kern2 .* ...
                    ((dev1)/bw_F).^p  .* ((dev2)/bw_F).^q);
                save_Q(nu_i,pq,tau_ii,tau_jj) = ss / nRawCovs(lag+1);
                if lag == 0
                    save_Q(nu_i,pq,tau_ii,tau_jj) = ss / ( censor.nGridTime * (bar_n2 - bar_n) );
                else
                    save_Q(nu_i,pq,tau_ii,tau_jj) = ss / ( (censor.nGridTime - lag) * bar_n^2 );
                end
            end

            % save S
            for pq = 1:6
                switch pq
                    case 1, p=0; q=0;
                    case 2, p=1; q=0;                                   
                    case 3, p=0; q=1;
                    case 4, p=1; q=1;
                    case 5, p=2; q=0;
                    case 6, p=0; q=2;    
                end
                ss = sum(...
                    kern1 .* ...
                    kern2 .* ...
                    ((dev1)/bw_F).^p  .* ((dev2)/bw_F).^q);
                if lag == 0
                    save_S(nu_i,pq,tau_ii,tau_jj) = ss / ( censor.nGridTime * (bar_n2 - bar_n) );
                else
                    save_S(nu_i,pq,tau_ii,tau_jj) = ss / ( (censor.nGridTime - lag) * bar_n^2 );
                end
            end
        end
    end            
end


%% complete save_Q and save_S for the negative lags
for lag=(-specDensity.q_bartlett+1):(-1)
    nu_i = lag + specDensity.q_bartlett; % the index of the real lag in for save_S and save_Q
    nu_i_symetric = -lag + specDensity.q_bartlett;

    for cv_jj = 1:censor.cv_K
        for cv_kk = 1:censor.cv_K

            % save Q
            for pq = 1:3
                switch pq % reverse p and q
                    case 1, qp = 1;
                    case 2, qp = 3;  
                    case 3, qp = 2;
                end
                save_Q(nu_i,pq,:,:) = squeeze(save_Q(nu_i_symetric,qp,:,:))';
            end

            % save S
            for pq = 1:6
                switch pq % reverse p and q
                    case 1, qp = 1;
                    case 2, qp = 3;
                    case 3, qp = 2;
                    case 4, qp = 4;
                    case 5, qp = 6;
                    case 6, qp = 5;
                end
                save_S(nu_i,pq,:,:) = squeeze(save_S(nu_i_symetric,qp,:,:))';
            end   
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%         PERFOM THE MANIPULATION TO EVALUATE THE DENSITY ON THE FREQ GRID             %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prepare S
S_00 = zeros(onb.nGridSmoother);
S_10 = zeros(onb.nGridSmoother);
S_01 = zeros(onb.nGridSmoother);
S_11 = zeros(onb.nGridSmoother);
S_20 = zeros(onb.nGridSmoother);
S_02 = zeros(onb.nGridSmoother);
for lag=(-specDensity.q_bartlett+1):(specDensity.q_bartlett-1)
    nu_i = lag + specDensity.q_bartlett;
    S_00 = S_00 + squeeze( save_S(nu_i,1,:,:) ) * (1-abs(lag)/specDensity.q_bartlett);
    S_10 = S_10 + squeeze( save_S(nu_i,2,:,:) ) * (1-abs(lag)/specDensity.q_bartlett);
    S_01 = S_01 + squeeze( save_S(nu_i,3,:,:) ) * (1-abs(lag)/specDensity.q_bartlett);
    S_11 = S_11 + squeeze( save_S(nu_i,4,:,:) ) * (1-abs(lag)/specDensity.q_bartlett);
    S_20 = S_20 + squeeze( save_S(nu_i,5,:,:) ) * (1-abs(lag)/specDensity.q_bartlett);
    S_02 = S_02 + squeeze( save_S(nu_i,6,:,:) ) * (1-abs(lag)/specDensity.q_bartlett);
end

% prepate some quantities
A_1 = S_20.*S_02 - S_11 .^2; % remember, all operations are pointwise
A_2 = S_10.*S_02 - S_01.*S_11;
A_3 = S_01.*S_20 - S_10.*S_11;
B = A_1.*S_00 - A_2.*S_10 - A_3.*S_01;

% estimate the spec density using the Bartlett's weights
specDensity.smoother = zeros(specDensity.nGridFreq, onb.nGridSmoother, onb.nGridSmoother); % this is in the grid for which I smoothed the lag covariances
for k = 1:specDensity.nGridFreq
    omega = specDensity.gridFreq(k); % real omega: -pi to pi       
    Q_00 = zeros(onb.nGridSmoother);
    Q_10 = zeros(onb.nGridSmoother);
    Q_01 = zeros(onb.nGridSmoother);
    for lag=(-specDensity.q_bartlett+1):(specDensity.q_bartlett-1)
        nu_i = lag + specDensity.q_bartlett;
        Q_00 = Q_00 + squeeze(save_Q(nu_i, 1, :, :)) * exp(-1i*lag*omega) * (1-abs(lag)/specDensity.q_bartlett);
        Q_10 = Q_10 + squeeze(save_Q(nu_i, 2, :, :)) * exp(-1i*lag*omega) * (1-abs(lag)/specDensity.q_bartlett);
        Q_01 = Q_01 + squeeze(save_Q(nu_i, 3, :, :)) * exp(-1i*lag*omega) * (1-abs(lag)/specDensity.q_bartlett);
    end

    % the spectral density
    m = 1/(2*pi) * specDensity.q_bartlett * (A_1.*Q_00 - A_2.*Q_10 - A_3.*Q_01) ./ B; % mind the fact that the sum of Batlett's weights is specDensity.q_bartlett
    specDensity.smoother(k,:,:) = (m+m')/2; % make sure it's self-adjoint
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the operator norms and traces-class norms
specDensity.norms = zeros(specDensity.nGridFreq,1);
specDensity.trace_class_norms = zeros(specDensity.nGridFreq,1);
for ii = 1:specDensity.nGridFreq
    
    % operator norm
    m = squeeze(specDensity.smoother(ii,:,:));
    specDensity.norms(ii) = norm( m );
    
    % trace class norm    
    specDensity.trace_class_norms(ii) = sum(abs(svd(m)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   positifiy the spectral density

% convert the spec density from onb.gridSmoother into the ONB
% positify the spectral densities, truncate the space corresponding to negative eigenvalues
specDensity.ONB = zeros(nGridFreq, onb.nBasis, onb.nBasis);
specDensity.ONB_positified = zeros(nGridFreq, onb.nBasis, onb.nBasis);
for k = 1:nGridFreq
    m = onb.smoothingMatrix_gridSmoother2ONB * squeeze(specDensity.smoother(k,:,:)) * onb.smoothingMatrix_gridSmoother2ONB';
    m = (m+m')/2; % make sure it's self-adjoint
    specDensity.ONB(k,:,:) = m;
    [V,D] = eig(m); % i.e. the decomposition:  m = V*D*V'
    diagD = real(diag(D)); % note that the harmonic eigenvalues must be always real
    diagD( diagD< 0 ) = 0; % truncate those elements that are less than (numerical) zero
    m = V*diag(diagD)*V';
    specDensity.ONB_positified(k,:,:) = (m+m')/2; % again, make sure it's self-adjoint  
end





end