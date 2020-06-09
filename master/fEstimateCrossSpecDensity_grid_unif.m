function crossSpecDensity = fEstimateCrossSpecDensity_grid_unif(x, data, mu_est, bw, nGridFreq, q_bartlett)

% define spectral weights
% weight_function = @(x) (abs(x)<=0.5) .* ( 1-6*x.^2 + 6*abs(x).^3 ) + (abs(x)>0.5) .* (2*(1-abs(x)).^3); % parzen weight
weight_function = @(x) 1-abs(x); % bartlett weight

% the mean of the TS "x"
x_mean = mean(x,'omitnan');

% define the data structure
crossSpecDensity = [];
crossSpecDensity.nGridFreq = nGridFreq;
crossSpecDensity.gridFreq = linspace(-pi,pi,crossSpecDensity.nGridFreq);

% shortcuts
onb = data.onb; % ONB
g = data.onb.gridSpace; %  for spatial grid

% selection rule for specDensity.q_bartlett
% mind the fact that there is also lag 0; so numOfLags=3 corresponds to lags 0,1,2
crossSpecDensity.q_bartlett = q_bartlett;
% specDensity.q_bartlett = floor( 2* data.nGridTime^(1/3));
% specDensity.q_bartlett = floor(data.nGridTime/2);


%%% calculate save_S save_Q

% here I'm going to save the tripple sum for Q's
% .. first index:  \nu = -specDensity.q_bartlett+1,...,+specDensity.q_bartlett-1
% .. second index: r:  1 = (r=0), 2 = (r=1)
% .. 3rd : evaluation on the onb.gridSmoother
save_Q = zeros(2*crossSpecDensity.q_bartlett-1, 2, onb.nGridSmoother);

% here I'm saving the S
% .. 1st:  \nu = -specDensity.q_bartlett+1,...,+specDensity.q_bartlett-1
% .. 2nd index: r:  1 = (r=0), 2 = (r=1), 3 = (r=2)
% .. 3rd : evaluation on the onb.gridSmoother
save_S = zeros(2*crossSpecDensity.q_bartlett-1, 3, onb.nGridSmoother);




for lag_real =(-crossSpecDensity.q_bartlett+1):(crossSpecDensity.q_bartlett-1)
    %disp(['calculating save Q and save S ... ',num2str(lag_real)])
    lag_indx = lag_real + crossSpecDensity.q_bartlett; % the index of the real lag in for save_S and save_Q
    
    % prepare the raw covariances for lag "lag"
    raws_G = zeros(data.n_maturities,1);
    raws_n = zeros(data.n_maturities,1);
    for m1 = 1:data.n_maturities
        
        
        % indeces of the maturities on teh spatial grid
        tau_1_indx = data.maturities_gs_indx(m1);
        
        % the vector of raw covariances, possibly with NaN elements
        shifted_start = max( lag_real+1, 1);
        shifted_end   = min( data.nGridTime + lag_real, data.nGridTime );
        x_start = max( 1, 1-lag_real );
        x_end   = min( data.nGridTime-lag_real, data.nGridTime );
                
        v = (data.yields( shifted_start:shifted_end ,m1) - mu_est.inSpace(tau_1_indx)) .* ...
            (x( x_start:x_end ) - x_mean );
        
        % calculate the mean of nonNaN elements and their number
        raws_G(m1) = mean(v, 'omitnan');
        raws_n(m1) = sum(~isnan(v));
        
    end
    
    % calculating Q^{(h)} and S^{(h)}
    for tau_ii = 1:onb.nGridSmoother
        tau_ii_inSpace = onb.gridSmoother(tau_ii); % the real position of tau_jj on [0,1] for the smoother
        
        % calculating Q^{(h)}(tau_ii_inSpace,tau_jj_inSpace) and S^{(h)}(tau_ii_inSpace,tau_jj_inSpace)
        
        % deviation from the maturities
        deviations = tau_ii_inSpace - g(data.maturities_gs_indx);
        
        % if there is at least one maturity close in both dimensions
        if (min(abs(deviations)) < bw)
            
            % spatial weights = Gaussian kernel
            W = normpdf(deviations/bw);
            
            % multiply by the weights corresponding to the number of available measurements
            W = W .* sum(raws_n);
            
            % save Q
            for pq = 1:3
                switch pq
                    case 1, r=0;
                    case 2, r=1;
                end
                save_Q(lag_indx,pq,tau_ii) = sum(...
                    raws_G' .*...
                    W .*...
                    (deviations/bw).^r );
            end
            
            % save S
            for pq = 1:6
                switch pq
                    case 1, r=0;
                    case 2, r=1;
                    case 3, r=2;
                end
                save_S(lag_indx,pq,tau_ii) = sum(...
                    W .*...
                    (deviations/bw).^r );
            end
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
for lag_real=(-crossSpecDensity.q_bartlett+1):(crossSpecDensity.q_bartlett-1)
    lag_indx = lag_real + crossSpecDensity.q_bartlett;
    S_0 = S_0 + squeeze( save_S(lag_indx,1,:,:) ) * (1-abs(lag_real)/crossSpecDensity.q_bartlett);
    S_1 = S_1 + squeeze( save_S(lag_indx,2,:,:) ) * (1-abs(lag_real)/crossSpecDensity.q_bartlett);
    S_2 = S_2 + squeeze( save_S(lag_indx,3,:,:) ) * (1-abs(lag_real)/crossSpecDensity.q_bartlett);
end


% estimate the spec density using the Bartlett's weights
crossSpecDensity.smoother = zeros(crossSpecDensity.nGridFreq, onb.nGridSmoother); % this is in the grid for which I smoothed the lag covariances
crossSpecDensity.smoother = zeros(crossSpecDensity.nGridFreq, onb.nBasis);
for k = 1:crossSpecDensity.nGridFreq
    omega = crossSpecDensity.gridFreq(k); % real omega: -pi to pi       
    Q_0 = zeros(onb.nGridSmoother,1);
    Q_1 = zeros(onb.nGridSmoother,1);
    for lag=(-crossSpecDensity.q_bartlett+1):(crossSpecDensity.q_bartlett-1)
        lag_index = lag + crossSpecDensity.q_bartlett;
        Q_0 = Q_0 + squeeze(save_Q(lag_index, 1, :)) * exp(-1i*lag*omega) * weight_function(abs(lag)/crossSpecDensity.q_bartlett);
        Q_1 = Q_1 + squeeze(save_Q(lag_index, 2, :)) * exp(-1i*lag*omega) * weight_function(abs(lag)/crossSpecDensity.q_bartlett);
    end

    % the spectral CROSS-density
    m = 1/(2*pi) * crossSpecDensity.q_bartlett * (Q_0.*S_2 - S_1.*Q_1) ./ (S_0.*S_2-S_1.^2); % mind the fact that the sum of Batlett's weights is specDensity.q_bartlett    
    crossSpecDensity.smoother(k,:) = m;
    crossSpecDensity.ONB(k,:) = onb.smoothingMatrix_gridSmoother2ONB * m;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculate the operator norms and traces-class norms
crossSpecDensity.norms = zeros(crossSpecDensity.nGridFreq,1);
crossSpecDensity.trace_class_norms = zeros(crossSpecDensity.nGridFreq,1);
for ii = 1:crossSpecDensity.nGridFreq
    
    % operator norm
    m = squeeze(crossSpecDensity.smoother(ii,:,:));
    crossSpecDensity.norms(ii) = norm( m );
    
    % trace class norm
    crossSpecDensity.trace_class_norms(ii) = sum(abs(svd(m)));
end






end