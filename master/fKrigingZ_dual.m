function response_est = fKrigingZ_dual( data_1, data_2, transfer_1, transfer_2, mu_ONB_1, mu_ONB_2, response_mean, dataFull_padding)

% prepara empty vector for output
nGridTime = size(data_1.est, 1) - 2*dataFull_padding;
response_est = response_mean * ones(1,nGridTime); % start with the response mean

% fill in the response_est vector
for t = 1:nGridTime
%     contributions = [];
    for lag_index = 1:(2*transfer_1.numOfLags+1)
        lag_real = lag_index - transfer_1.numOfLags-1;
        
        if (t-lag_real+dataFull_padding >= 1) && (t-lag_real+dataFull_padding <= size(data_1.est, 1))
            contribution_1 = transfer_1.ONB(lag_index,:) * ( data_1.est( t-lag_real+dataFull_padding ,:)' - mu_ONB_1 );            
            contribution_2 = transfer_2.ONB(lag_index,:) * ( data_2.est( t-lag_real+dataFull_padding ,:)' - mu_ONB_2 );            
%             contributions = [contributions contribution];
            response_est(t) = response_est(t) + contribution_1 + contribution_2;        
        end
    end    
end

end