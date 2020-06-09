function response_est = fKrigingZ(krigingX_padded,transfer,dataFull_padding)

% prepara empty vector for output
nGridTime = size(krigingX_padded.est, 1) - 2*dataFull_padding;
response_est = zeros(1,nGridTime);

% fill in the response_est vector
for t = 1:nGridTime
    for lag_index = 1:(2*transfer.numOfLags+1)
        lag_real = lag_index - transfer.numOfLags-1;
        response_est(t) = response_est(t) + transfer.ONB(lag_index,:) * krigingX_padded.est( t-lag_real+dataFull_padding ,:)';
    end    
end

end