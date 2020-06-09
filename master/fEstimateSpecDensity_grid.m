function specDensity = fEstimateSpecDensity_grid(data,mu_est,bw_1,bw_2,nGridFreq, q_bartlett)
% this version has different rule for Bartletts span parameter

specDensity = [];
specDensity.nGridFreq = nGridFreq;
specDensity.gridFreq = linspace(-pi,pi,specDensity.nGridFreq);

% shortcuts
onb = data.onb; % ONB
g = data.onb.gridSpace; %  for spatial grid

% selection rule for specDensity.q_bartlett
% mind the fact that there is also lag 0; so numOfLags=3 corresponds to lags 0,1,2
specDensity.q_bartlett = q_bartlett;
% specDensity.q_bartlett = floor( 2* data.nGridTime^(1/3));
% specDensity.q_bartlett = floor(data.nGridTime/2);

%%% estimate the average number of measurements per curve
bar_n = mean( sum( ~isnan(data.yields), 2 ));
bar_n2 = mean( sum( ~isnan(data.yields), 2 ).^2);

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
    
    % prepare the raw covariances for lag "lag"
    raws_G = zeros(data.n_maturities);
    raws_n = zeros(data.n_maturities);
    for m1 = 1:data.n_maturities
        for m2 = 1:data.n_maturities
            
            % discard the diagonal if lag==0
            if (lag ~= 0) || (m1 ~= m2)

                % indeces of the maturities on teh spatial grid
                tau_1_indx = data.maturities_gs_indx(m1);
                tau_2_indx = data.maturities_gs_indx(m2);

                % the vector of raw covariances, possibly with NaN elements
                v = (data.yields((lag+1):end,m1) - mu_est.inSpace(tau_1_indx)) .* (data.yields(1:(data.nGridTime-lag),m2) - mu_est.inSpace(tau_2_indx));
                
                % calculate the mean of nonNaN elements and their number
                raws_G(m1,m2) = mean(v, 'omitnan');
                raws_n(m1,m2) = sum(~isnan(v));
            end
        end
    end

    % calculating Q^{(h)} and S^{(h)}
    for tau_ii = 1:onb.nGridSmoother
        tau_ii_inSpace = onb.gridSmoother(tau_ii); % the real position of tau_jj on [0,1] for the smoother
        for tau_jj = 1:onb.nGridSmoother                            
            tau_jj_inSpace = onb.gridSmoother(tau_jj); % the real position of tau_kk on [0,1] for the smoother

            % calculating Q^{(h)}(tau_ii_inSpace,tau_jj_inSpace) and S^{(h)}(tau_ii_inSpace,tau_jj_inSpace)
            
            % deviation from the maturities
            deviations_1 = tau_ii_inSpace - g(data.maturities_gs_indx);
            deviations_2 = tau_jj_inSpace - g(data.maturities_gs_indx);
            deviation_matrix_1 = repmat(deviations_1,  length(data.maturities_gs_indx), 1);
            deviation_matrix_2 = repmat(deviations_2', 1, length(data.maturities_gs_indx));

            % if there is at least one maturity close in both dimensions
            if (min(abs(deviations_1)) < bw_1(tau_ii,tau_jj)) && (min(abs(deviations_2)) < bw_2(tau_ii,tau_jj))
                
                % spatial weights = Gaussian kernel
                W = normpdf(deviation_matrix_1(:)/bw_1(tau_ii,tau_jj)) .* normpdf(deviation_matrix_2(:)/bw_2(tau_ii,tau_jj));
                
                % multiply by the weights corresponding to the number of available measurements
                W = W .* sum(raws_n(:));
                                
                % save Q
                for pq = 1:3
                    switch pq
                        case 1, p=0; q=0;
                        case 2, p=1; q=0;
                        case 3, p=0; q=1;
                    end
                    ss = sum(...
                        raws_G(:) .*...
                        W .*...
                        (deviation_matrix_1(:)/bw_1(tau_ii,tau_jj)).^p .*...
                        (deviation_matrix_2(:)/bw_2(tau_ii,tau_jj)).^q );
                    if lag == 0
                        save_Q(nu_i,pq,tau_ii,tau_jj) = ss / ( data.nGridTime * (bar_n2 - bar_n) );
                    else
                        save_Q(nu_i,pq,tau_ii,tau_jj) = ss / ( (data.nGridTime - lag) * bar_n^2 );
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
                        W .*...
                        (deviation_matrix_1(:)/bw_1(tau_ii,tau_jj)).^p .*...
                        (deviation_matrix_2(:)/bw_2(tau_ii,tau_jj)).^q );
                    if lag == 0
                        save_S(nu_i,pq,tau_ii,tau_jj) = ss / ( data.nGridTime * (bar_n2 - bar_n) );
                    else
                        save_S(nu_i,pq,tau_ii,tau_jj) = ss / ( (data.nGridTime - lag) * bar_n^2 );
                    end
                    
                end                
            end            
        end
    end            
end


%% complete save_Q and save_S for the negative lags
for lag=(-specDensity.q_bartlett+1):(-1)
    nu_i = lag + specDensity.q_bartlett; % the index of the real lag in for save_S and save_Q
    nu_i_symetric = -lag + specDensity.q_bartlett;
    
    
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