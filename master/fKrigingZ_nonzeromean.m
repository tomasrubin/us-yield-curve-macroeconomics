function response_est = fKrigingZ_nonzeromean(krigingX_padded,transfer, mu_ONB, response_mean ,dataFull_padding)

% prepara empty vector for output
nGridTime = size(krigingX_padded.est, 1) - 2*dataFull_padding;
response_est = response_mean * ones(1,nGridTime); % start with the response mean

% fill in the response_est vector
for t = 1:nGridTime
%     contributions = [];
    for lag_index = 1:(2*transfer.numOfLags+1)
        lag_real = lag_index - transfer.numOfLags-1;
        
        if (t-lag_real+dataFull_padding >= 1) && (t-lag_real+dataFull_padding <= size(krigingX_padded.est, 1))
            contribution = transfer.ONB(lag_index,:) * ( krigingX_padded.est( t-lag_real+dataFull_padding ,:)' - mu_ONB );            
%             contributions = [contributions contribution];
            response_est(t) = response_est(t) + contribution;        
        end
    end    
end

end