function [response, response_wo_noise, transfer_true, specCrossdensity_true] = fSimulate_response_ONB_v3(onb,response_process,response_process_shape,dataFullONB_padded,dataFull_padding,specDensity_true)
% inputs
%   onb ... struct for ONB
%   response_process ... integer - what response equation I want to use
%   pars_response_process ... matrix(nFunctions,nGridSpace) of pre-generated random functions
%   dataFull_padded .... matrix - data including padding before and after the true data sample
%   dataFull_padding ... integer - how much has been padded in the matrix above
%   specCrossdensity_true ... struct - required for calculating the cross-spectral density

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prepare some constants


% all response equations are up to 10 lags each sides, thus 21 parameters in total
transfer_true = [];
transfer_true.numOfLags = 10;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        set up the transfer parameters

g = onb.gridSpace;
switch response_process_shape
    case 1
        disp("fSimulate_response_ONB_v3: pars_response_process = 1")    
        pars_response_process = 1;
    case 2
        disp("fSimulate_response_ONB_v3: pars_response_process =  cos( 2*g*2*pi )")
        pars_response_process = cos( 2*g*2*pi );
    case 3
        disp("fSimulate_response_ONB_v3: pars_response_process =  sin( g*2*pi )")
        pars_response_process = sin( g*2*pi );
    case 4
        disp("fSimulate_response_ONB_v3: pars_response_process =  ( -g+0.5).* exp( -20*( -g+0.5).^2 )")
        pars_response_process = ( -g+0.5).* exp( -20*( -g+0.5).^2 );
    
    otherwise % unknown process => stop the script
        fprintf('Error: unknown response_process_shape\n')
        return;
end

transfer_true.inSpace = zeros(transfer_true.numOfLags*2+1, onb.nGridSpace);
switch response_process
    case 1
        transfer_true.inSpace(transfer_true.numOfLags+1,:) = pars_response_process; % b_0
        transfer_true.inSpace(transfer_true.numOfLags+2,:) = pars_response_process; % b_1
        transfer_true.sigma = sqrt(0.001);
        
    case 2
        transfer_true.inSpace(transfer_true.numOfLags+1,:) = pars_response_process; % b_0
        transfer_true.inSpace(transfer_true.numOfLags+4,:) = pars_response_process; % b_3
        transfer_true.sigma = sqrt(0.001);
        
    case 3
        transfer_true.inSpace(transfer_true.numOfLags+1,:) =   1 * pars_response_process; % b_0
        transfer_true.inSpace(transfer_true.numOfLags+2,:) = 0.9 * pars_response_process; % b_1
        transfer_true.inSpace(transfer_true.numOfLags+3,:) = 0.7 * pars_response_process; % b_2
        transfer_true.inSpace(transfer_true.numOfLags+4,:) = 0.5 * pars_response_process; % b_3
        transfer_true.inSpace(transfer_true.numOfLags+5,:) = 0.3 * pars_response_process; % b_4
        transfer_true.inSpace(transfer_true.numOfLags+6,:) = 0.1 * pars_response_process; % b_5
        transfer_true.sigma = sqrt(0.001);
%         
%     case 11
%         transfer_true.inSpace(transfer_true.numOfLags+1,:) = pars_response_process(1,:); % b_0
%         transfer_true.inSpace(transfer_true.numOfLags+2,:) = pars_response_process(2,:); % b_1
%         transfer_true.sigma = sqrt(0.5);
%         
%     case 12
%         transfer_true.inSpace(transfer_true.numOfLags+1,:) = pars_response_process(1,:); % b_0
%         transfer_true.inSpace(transfer_true.numOfLags+4,:) = pars_response_process(2,:); % b_3
%         transfer_true.sigma = sqrt(0.5);
%         
%     case 13
%         transfer_true.inSpace(transfer_true.numOfLags+1,:) =   1 * pars_response_process(1,:); % b_0
%         transfer_true.inSpace(transfer_true.numOfLags+2,:) = 0.9 * pars_response_process(2,:); % b_1
%         transfer_true.inSpace(transfer_true.numOfLags+3,:) = 0.7 * pars_response_process(3,:); % b_2
%         transfer_true.inSpace(transfer_true.numOfLags+4,:) = 0.5 * pars_response_process(4,:); % b_3
%         transfer_true.inSpace(transfer_true.numOfLags+5,:) = 0.3 * pars_response_process(5,:); % b_4
%         transfer_true.inSpace(transfer_true.numOfLags+6,:) = 0.1 * pars_response_process(6,:); % b_5
%         transfer_true.sigma = sqrt(0.5);
        
        
    otherwise % unknown process => stop the script
        fprintf('Error: unknown process\n')
        return;
end

% convert all into ONB
transfer_true.ONB = nan(transfer_true.numOfLags*2+1, onb.nBasis);
for lag=1:(transfer_true.numOfLags*2+1)
    transfer_true.ONB(lag,:) = onb.onbMatrix' * transfer_true.inSpace(lag,:)' / onb.nGridSpace;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        calculate the true spectral response operator

% prepera the spectral data structure
transfer_true.nGridFreq = specDensity_true.nGridFreq;
transfer_true.gridFreq = linspace(-pi,pi, transfer_true.nGridFreq);
transfer_true.specTransfer = zeros(transfer_true.nGridFreq, onb.nBasis);

% sum up for each frequency
for k = 1:transfer_true.nGridFreq
    omega = transfer_true.gridFreq(k);
    
    % sum up over all lags
    m = zeros(1,onb.nBasis); % prepare zeros matrix
    for lag = 1:(transfer_true.numOfLags*2+1)
        lag_real = lag-transfer_true.numOfLags-1; % real lag, i.e. -numOfLags,...,-1,0,1,...,numOfLags
        m = m + transfer_true.ONB(lag,:) * exp(-1i*lag_real*omega);
    end    
    transfer_true.specTransfer(k,:) = m;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        calculate the cross-spectral density
% formula: F_\omega^{YX} = B_\omega o F_\omega^{X}

% prepare the data structure
specCrossdensity_true = [];
specCrossdensity_true.nGridFreq = transfer_true.nGridFreq;
specCrossdensity_true.gridFreq = linspace(-pi,pi, specCrossdensity_true.nGridFreq);
specCrossdensity_true.ONB = zeros(specCrossdensity_true.nGridFreq, onb.nBasis);

% cycle across all frequencies
for k = 1:specCrossdensity_true.nGridFreq
    specCrossdensity_true.ONB(k,:) = transfer_true.specTransfer(k,:) * squeeze(specDensity_true.ONB(k,:,:));
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        calculate the response process Z

% extract the size of the dataset
[nGridTime_padded, ~] = size(dataFullONB_padded);
nGridTime = nGridTime_padded - 2*dataFull_padding;

% output vector
response_wo_noise = zeros(nGridTime,1);

% cycle over all time points of real dataset (without padding)
for t = 1:nGridTime
    
    % contributions from the main FTS X
    for lag = 1:(transfer_true.numOfLags*2+1)
        lag_real = lag-transfer_true.numOfLags-1;
        response_wo_noise(t) = response_wo_noise(t) + dataFullONB_padded(t+dataFull_padding-lag_real,:) * transfer_true.ONB(lag,:)' ;
    end    
    
    
end

% model noise
response = response_wo_noise + transfer_true.sigma*randn(size(response_wo_noise));


end