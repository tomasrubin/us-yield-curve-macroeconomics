function cov0_est = fSmootherCov0_grid(data, bw_adjustment, mu_est)

% the output variables
cov0_est = [];
cov0_est.gridSmoother = nan(data.onb.nGridSmoother);

% create raw covariances
raws_G = zeros(data.n_maturities);
raws_n = zeros(data.n_maturities);
for m1 = 1:data.n_maturities
    for m2 = 1:data.n_maturities
        
        if (m1 ~= m2) % discard the diagonal

            % indeces of the maturities on teh spatial grid
            tau_1_indx = data.maturities_gs_indx(m1);
            tau_2_indx = data.maturities_gs_indx(m2);

            % create raw covariances
            v = (data.yields(:,m1) - mu_est.inSpace(tau_1_indx)) .* (data.yields(:,m2) - mu_est.inSpace(tau_2_indx));
            raws_G(m1,m2) = mean(v, 'omitnan');
            raws_n(m1,m2) = sum(~isnan(v));

        end
    end
end

% distance to closest maturity determines the BW
distance = nan(1,data.onb.nGridSmoother);
distance_old = 0.01; % the minimal BW from which I start
for ii = 1:data.onb.nGridSmoother
    tau_indx = data.onb.gs_indx(ii);
    distance_new = min(abs(data.onb.gridSpace(tau_indx) - data.onb.gridSpace(data.maturities_gs_indx)));
    if distance_new < distance_old
        distance_new = distance_old; % if the new distance is less => use the old distance
    else
        distance_old = distance_new; % if the new distance is greater => use the new one, and save it for next time
    end
    distance(ii) = distance_new;    
end

% determine the bandwidths
cov0_est.bw_2 = repmat(distance,data.onb.nGridSmoother,1) * bw_adjustment;
cov0_est.bw_1 = repmat(distance',1,data.onb.nGridSmoother) * bw_adjustment;


% perform the smoother
for ii = 1:data.onb.nGridSmoother
    for jj = 1:data.onb.nGridSmoother
        tau1_indx = data.onb.gs_indx(ii);
        tau2_indx = data.onb.gs_indx(jj);
        cov0_est.gridSmoother(ii,jj) = fSmootherSurface_grid_point( tau1_indx, tau2_indx, raws_G, raws_n, data.maturities_gs_indx,...
            cov0_est.bw_1(ii,jj), cov0_est.bw_2(ii,jj), data.onb);
    end
end


% convert to the basis
cov0_est.ONB = data.onb.smoothingMatrix_gridSmoother2ONB * cov0_est.gridSmoother * data.onb.smoothingMatrix_gridSmoother2ONB';
cov0_est.inSpace = data.onb.onbMatrix * cov0_est.ONB * data.onb.onbMatrix';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%   function to calculate the CV score for a given candidate bandwidth
%     function mse = sum_of_squares(m,bw)
%
%         mse = 0;
%         % for each K partition
%         for k = 1:data.n_cv_k
%
%             % fit the smoother on all except "k" and except maturity "m"
%             Y_train = data.yields;
%             Y_train(:,m) = NaN; % except maturity "m"
%             Y_train( find(data.cv_k == k) ,:) = NaN; % except "k"
%
%             % perform the smoother
%             mu_est = fSmootherLine_grid_point( data.maturities_gs_indx(m), Y_train, data.maturities_gs_indx, bw, data.onb);
%
%             % compare with the test partition
%             Y_test = data.yields( data.cv_k == k, m ); % only test on partition "k"
%
%             if ~isempty(Y_test)
%                 mse = mse + 1/data.n_cv_k * mean((mu_est - Y_test).^2 ,'omitnan');
%             end
%         end
%
%
%     end

end