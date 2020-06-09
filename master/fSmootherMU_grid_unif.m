function mu_est = fSmootherMU_grid_unif(data, bw)

% the output variables
mu_est = [];
mu_est.gridSmoother = nan(data.onb.nGridSmoother,1);
mu_est.bw = ones(data.onb.nGridSmoother,1) * bw;

% % distance to closest maturity determines the BW
% distance = nan(1,data.onb.nGridSmoother);
% distance_old = 0.01; % the minimal BW from which I start
% for ii = 1:data.onb.nGridSmoother
%     tau_indx = data.onb.gs_indx(ii);
%     distance_new = min(abs(data.onb.gridSpace(tau_indx) - data.onb.gridSpace(data.maturities_gs_indx)));
%     if distance_new < distance_old
%         distance_new = distance_old; % if the new distance is less => use the old distance
%     else
%         distance_old = distance_new; % if the new distance is greater => use the new one, and save it for next time
%     end
%     distance(ii) = distance_new;
%     mu_est.bw(ii) = distance(ii) * bw_adjustment;
% end

% perform the smoother
for ii = 1:data.onb.nGridSmoother
    tau_indx = data.onb.gs_indx(ii);    
%     bw = (1-data.onb.gridSpace(tau_indx)) * bw_start + data.onb.gridSpace(tau_indx) * bw_end;    
%     mu_est.bw(ii) = bw(ii);
    mu_est.gridSmoother(ii) = fSmootherLine_grid_point( tau_indx, data.yields, data.maturities_gs_indx, mu_est.bw(ii), data.onb);
end


% convert to the basis
mu_est.ONB = data.onb.smoothingMatrix_gridSmoother2ONB * mu_est.gridSmoother;
mu_est.inSpace = data.onb.onbMatrix * mu_est.ONB;


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