function [transfer_1,transfer_2] = fEstimateTransfer_truncation_dual(onb,specDensity_1,specDensity_2,specDensity_1_2,specCrossDensity_1,specCrossDensity_2, truncationParameter_1, truncationParameter_2)

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
    
    % spectral truncation 1
    M_1 = (squeeze(specDensity_1.ONB(k,:,:)) + squeeze(specDensity_1.ONB(k,:,:))')/2;    
    [V_1,D_1] = eig(M_1);
    diagD_1 = diag(D_1);
    truncation_number_eigen_1 = sum(diagD_1 > truncationParameter_1); % count how many eigenvalues are above the threshold
    transfer_1.truncationEigen(k) = truncation_number_eigen_1;
    [eigenvalues_1,indeces_ordering_1] = sort(diagD_1, 'descend');
    eigenvalues_1 = eigenvalues_1(1:truncation_number_eigen_1); % cut of eigenvalues beyond the truncation level
    V_1 = V_1(:,indeces_ordering_1); % reorder the eigenfunctions
    V_1 = V_1(:,1:truncation_number_eigen_1); % take only those eigenfunction corresponding to the correct truncation
    
    % spectral truncation 2
    M_2 = (squeeze(specDensity_2.ONB(k,:,:)) + squeeze(specDensity_2.ONB(k,:,:))')/2;    
    [V_2,D_2] = eig(M_2);
    diagD_2 = diag(D_2);
    truncation_number_eigen_2 = sum(diagD_2 > truncationParameter_2); % count how many eigenvalues are above the threshold
    transfer_2.truncationEigen(k) = truncation_number_eigen_2;
    [eigenvalues_2,indeces_ordering_2] = sort(diagD_2, 'descend');
    eigenvalues_2 = eigenvalues_2(1:truncation_number_eigen_2); % cut of eigenvalues beyond the truncation level
    V_2 = V_2(:,indeces_ordering_2); % reorder the eigenfunctions
    V_2 = V_2(:,1:truncation_number_eigen_2); % take only those eigenfunction corresponding to the correct truncation
    
    % join basis for V_1, V_2
    V_12 = zeros(onb.nBasis*2, truncation_number_eigen_1 + truncation_number_eigen_2);
    V_12(1:onb.nBasis, 1:truncation_number_eigen_1) = V_1;
    V_12((onb.nBasis+1):end, (truncation_number_eigen_1+1):end) = V_2;
    
    % express the joint spectral density in the reduced basis
    M = nan(truncation_number_eigen_1 + truncation_number_eigen_2);
    
    % blocks corresponding to the marginal spectral densities
    M( 1:truncation_number_eigen_1,1:truncation_number_eigen_1 ) = diag(eigenvalues_1);
    M( (truncation_number_eigen_1+1):end,(truncation_number_eigen_1+1):end ) = diag(eigenvalues_2);
    
    % cross density between 1 and 2
    M_12 = squeeze(specDensity_1_2.ONB(k,:,:));    
    M_12_in_reduced_basis = V_1' * M_12 * V_2;
    M( 1:truncation_number_eigen_1, (truncation_number_eigen_1+1):end ) = M_12_in_reduced_basis;
    M( (truncation_number_eigen_1+1):end, 1:truncation_number_eigen_1 ) = M_12_in_reduced_basis';
    
    % inversion in reduced basis
    M_inverted = V_12 * inv(M) * V_12';
    
    % pre-multiply by the spectral crossdensity
    x = [specCrossDensity_1.ONB(k,:), specCrossDensity_2.ONB(k,:)] * M_inverted;
    
    % extract the spectral transfer function for both components
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