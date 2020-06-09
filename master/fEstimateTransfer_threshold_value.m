function transfer = fEstimateTransfer_threshold_value(onb,specDensity,specCrossDensity, threshold_value)

% prepere the spectral data structure
transfer.nGridFreq = specDensity.nGridFreq;
transfer.gridFreq = linspace(-pi,pi, transfer.nGridFreq);
transfer.specTransfer = zeros(transfer.nGridFreq, onb.nBasis);
transfer.truncationEigen = nan(transfer.nGridFreq,1);

%% spectral transfer function
for k=1:transfer.nGridFreq
        
    % spectral truncation of M
    M = (squeeze(specDensity.ONB(k,:,:)) + squeeze(specDensity.ONB(k,:,:))')/2;
    [V,D] = eig(M);
    diagD = diag(D);
    truncationEigen = sum(diagD > threshold_value); % count how many eigenvalues are above the threshold
    transfer.truncationEigen(k) = truncationEigen;
    [B,I] = sort(diagD, 'descend');
    B(1:truncationEigen) = 1./B(1:truncationEigen); % inversion
    B((truncationEigen+1):end) = 0;
    
    % get the spectral transfer function
    
    transfer.specTransfer(k,:) = specCrossDensity.ONB(k,:) * V(:,I)*diag(B)*V(:,I)';
end

%% transfer functional --- integrate out the spectral transfer function
% all response equations are up to 10 lags each sides, thus 21 parameters in total
transfer.numOfLags = 30;

transfer.ONB = zeros(2*transfer.numOfLags+1,onb.nBasis);
for lag = 1:(2*transfer.numOfLags+1)
    lag_real = lag-transfer.numOfLags-1;
    
    for k=1:transfer.nGridFreq
        omega=transfer.gridFreq(k);
        transfer.ONB(lag,:) = transfer.ONB(lag,:) + 1/ transfer.nGridFreq * transfer.specTransfer(k,:) * exp(1i*lag_real*omega);
    end
end

% delete numerically non-zero imaginary part
transfer.ONB = real(transfer.ONB);

end