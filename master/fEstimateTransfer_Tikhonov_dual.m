function [transfer_1,transfer_2] = fEstimateTransfer_Tikhonov_dual(onb,specDensity_1,specDensity_2,specDensity_1_2,specCrossDensity_1,specCrossDensity_2, rho_1, rho_2)

% prepere the spectral data structure
transfer_1 = [];
transfer_1.nGridFreq = specDensity_1.nGridFreq;
transfer_1.gridFreq = linspace(-pi,pi, transfer_1.nGridFreq);
transfer_1.specTransfer = zeros(transfer_1.nGridFreq, onb.nBasis);
transfer_2 = [];
transfer_2.nGridFreq = specDensity_1.nGridFreq;
transfer_2.gridFreq = linspace(-pi,pi, transfer_1.nGridFreq);
transfer_2.specTransfer = zeros(transfer_1.nGridFreq, onb.nBasis);

%% spectral transfer function
for k=1:transfer_1.nGridFreq
    
    % spectral regularization of M
    M = zeros( onb.nBasis*2 );
    M(1:onb.nBasis,1:onb.nBasis) = squeeze(specDensity_1.ONB(k,:,:)) + rho_1 * eye(onb.nBasis);
    M(1:onb.nBasis,(onb.nBasis+1):end) = squeeze(specDensity_1_2.ONB(k,:,:))';
    M((onb.nBasis+1):end,1:onb.nBasis) = squeeze(specDensity_1_2.ONB(k,:,:));
    M((onb.nBasis+1):end,(onb.nBasis+1):end) = squeeze(specDensity_2.ONB(k,:,:)) + rho_2 * eye(onb.nBasis);
    
    M = (M+M')/2;
    
	x = [specCrossDensity_1.ONB(k,:), specCrossDensity_2.ONB(k,:)]  / M;
%     x = [specCrossDensity_1.ONB(k,:), zeros(size(specCrossDensity_2.ONB(k,:)))]  / M;
    transfer_1.specTransfer(k,:) = x(1:onb.nBasis);
    transfer_2.specTransfer(k,:) = x((onb.nBasis+1):end);
end



%% filter coefficients --- integrate out the spectral transfer function
% all response equations are up to 10 lags each sides, thus 21 parameters in total
transfer_1.numOfLags = 10;
transfer_2.numOfLags = transfer_1.numOfLags;

transfer_1.ONB = zeros(2*transfer_1.numOfLags+1,onb.nBasis);
transfer_2.ONB = zeros(2*transfer_1.numOfLags+1,onb.nBasis);
for lag = 1:(2*transfer_1.numOfLags+1)
    lag_real = lag-transfer_1.numOfLags-1;
    
    for k=1:transfer_1.nGridFreq
        omega=transfer_1.gridFreq(k);
        transfer_1.ONB(lag,:) = transfer_1.ONB(lag,:) + 1/ transfer_1.nGridFreq * transfer_1.specTransfer(k,:) * exp(1i*lag_real*omega);
        transfer_2.ONB(lag,:) = transfer_2.ONB(lag,:) + 1/ transfer_2.nGridFreq * transfer_2.specTransfer(k,:) * exp(1i*lag_real*omega);
    end
end

% delete numerically non-zero imaginary part
transfer_1.ONB = real(transfer_1.ONB);
transfer_2.ONB = real(transfer_2.ONB);

end