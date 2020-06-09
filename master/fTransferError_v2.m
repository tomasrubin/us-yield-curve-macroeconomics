function [rmse, mse] = fTransferError_v2(transfer,transfer_true)

nBasis = length(transfer_true.ONB(1,:));

% identify lags with non-zero
nonzero_transfers = zeros(2*transfer_true.numOfLags+1,1);
for lag_index = 1:(2*transfer_true.numOfLags+1)
    if max(abs(transfer_true.inSpace(lag_index,:))) > 1e-6 % if it's numerically nonzero, save it as a nonzero transfer functional
        nonzero_transfers(lag_index) = 1;
    end
end

% prepare the structure for RMSE of nonzero functionals
numOfFunctionals = sum(nonzero_transfers);
rmse_transfers = zeros(numOfFunctionals,1);

% cycle over all nonzero functionals
iFunctional = 1;
for lag = -transfer_true.numOfLags:transfer_true.numOfLags    
    if nonzero_transfers(lag + transfer_true.numOfLags + 1)
        rmse_transfers(iFunctional) = ...
            norm( transfer.ONB( lag + transfer.numOfLags + 1 ,:) -...
            transfer_true.ONB(lag + transfer_true.numOfLags + 1,:) )...
            / norm(transfer_true.ONB(lag + transfer_true.numOfLags + 1,:));
        iFunctional = iFunctional + 1;
    end    
end

rmse = mean(rmse_transfers);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% absolute error
mse = 0;
% cycle over all filters for which I have an estimate
for lag = -transfer.numOfLags:transfer.numOfLags
    
    truth_ONB = zeros( 1, nBasis );
    if abs(lag) <= transfer_true.numOfLags
        truth_ONB = transfer_true.ONB( transfer_true.numOfLags + 1 + lag, : );
    end
    
    mse = mse + norm(transfer.ONB( transfer.numOfLags + 1 + lag, :) - truth_ONB)^2;
    
end

end