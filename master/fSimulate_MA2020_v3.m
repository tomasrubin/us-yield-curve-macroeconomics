function [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covs_true, mu_true, specDensity_true, MA_order, sigma] = fSimulate_MA2020_v3(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding)

% some constants
cv_K = 10; % the number of partitions for K-fold CV


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              define sigma and MA_order

switch simulated_process
           
    case 3 %% MA(2)        
        MA_order = 2;
        
    case 4 %% MA(4)        
        MA_order = 4;
        
    case 5 %% MA(8)        
        MA_order = 8;
    
    otherwise % unknown process => stop the script
        fprintf('Error: unknown process\n')
        return;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              define the parameters of the MA process
%% set up the parameters
g = onb.gridSpace;

% covariance of stochastic term
stochCov_inSpace = ...
      1 * sin(   g*2*pi ).*sin(   g'*2*pi ) +...
    0.6 * cos(   g*2*pi ).*cos(   g'*2*pi ) +...
    0.3 * sin( 2*g*2*pi ).*sin( 2*g'*2*pi ) +...
    0.1 * cos( 2*g*2*pi ).*cos( 2*g'*2*pi ) +...
    0.1 * sin( 3*g*2*pi ).*sin( 3*g'*2*pi ) +...
    0.1 * cos( 3*g*2*pi ).*cos( 3*g'*2*pi ) +...
    0.05* sin( 4*g*2*pi ).*sin( 4*g'*2*pi ) +...
    0.05* cos( 4*g*2*pi ).*cos( 4*g'*2*pi ) +...
    0.05* sin( 5*g*2*pi ).*sin( 5*g'*2*pi ) +...
    0.05* cos( 5*g*2*pi ).*cos( 5*g'*2*pi );
par_covE = onb.onbMatrix' * stochCov_inSpace * onb.onbMatrix / onb.nGridSpace^2;

% MA coeficients
par_A1_inSpace = sin((g'+g)*2*pi);
par_A2_inSpace = sin(((1-g)'+g)*2*pi);
par_A3_inSpace = sin(((1-g)'+(1-g))*2*pi);
par_A4_inSpace = sin((g'+(1-g))*2*pi);
par_A1 = onb.onbMatrix' * par_A1_inSpace * onb.onbMatrix / onb.nGridSpace^2;
par_A2 = onb.onbMatrix' * par_A2_inSpace * onb.onbMatrix / onb.nGridSpace^2;
par_A3 = onb.onbMatrix' * par_A3_inSpace * onb.onbMatrix / onb.nGridSpace^2;
par_A4 = onb.onbMatrix' * par_A4_inSpace * onb.onbMatrix / onb.nGridSpace^2;

% normalize them
par_A1 = par_A1 / norm(par_A1) * 0.8;
par_A2 = par_A2 / norm(par_A2) * 0.6;
par_A3 = par_A3 / norm(par_A3) * 0.4;
par_A4 = par_A4 / norm(par_A4) * 0.2;

% mean function
% mu_true_inSpace = 4 * sin(onb.gridSpace*pi*1.5);
mu_true_inSpace = zeros(size(onb.gridSpace));
mu_true = onb.onbMatrix' * mu_true_inSpace' / onb.nGridSpace;





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              true auto-covariances

% pre-alocate space
covs_true = [];

% save how many lags I have and pre-alocate space
covs_true.numOfLags = nGridTime - 1;
covs_true.ONB = zeros(covs_true.numOfLags+1, onb.nBasis, onb.nBasis);
covs_true.ONB_norms = zeros(covs_true.numOfLags+1,1);


% save how many lags I have and pre-alocate space
covs_true.numOfLags = nGridTime - 1;
covs_true.ONB = zeros(covs_true.numOfLags+1, onb.nBasis, onb.nBasis);
covs_true.ONB_norms = zeros(covs_true.numOfLags+1,1);

if (MA_order == 2)
    covs_true.ONB(1,:,:) = par_covE + ...
        par_A1*par_covE*par_A1' +...
        par_A2*par_covE*par_A2';
    covs_true.ONB(2,:,:) = par_A1*par_covE + par_A2*par_covE*par_A1';
    covs_true.ONB(3,:,:) = par_A2*par_covE;
end
if (MA_order == 4)
    covs_true.ONB(1,:,:) = par_covE +...
        par_A1*par_covE*par_A1' +...
        par_A2*par_covE*par_A2' +...
        par_A3*par_covE*par_A3' +...
        par_A4*par_covE*par_A4';
    covs_true.ONB(2,:,:) = ...
        par_A1*par_covE +...
        par_A2*par_covE*par_A1' + ...
        par_A3*par_covE*par_A2' + ...
        par_A4*par_covE*par_A3';
    covs_true.ONB(3,:,:) = ...
        par_A2*par_covE + ...
        par_A3*par_covE*par_A1' + ...
        par_A4*par_covE*par_A2';
    covs_true.ONB(4,:,:) = ...
        par_A3*par_covE +...
        par_A4*par_covE*par_A1';
    covs_true.ONB(5,:,:) = par_A4*par_covE;
end
if (MA_order == 8)
    
    covs_true.ONB(1,:,:) = par_covE + ...
        par_A1*par_covE*par_A1' + ...
        par_A2*par_covE*par_A2' +...
        par_A3*par_covE*par_A3' + ...
        par_A4*par_covE*par_A4' +...
        par_A1*par_covE*par_A1' + ...
        par_A2*par_covE*par_A2' +...
        par_A3*par_covE*par_A3' + ...
        par_A4*par_covE*par_A4';
    covs_true.ONB(2,:,:) = ...
        par_A1*par_covE + ...
        par_A2*par_covE*par_A1' +...
        par_A3*par_covE*par_A2' + ...
        par_A4*par_covE*par_A3' +...
        par_A1*par_covE*par_A4' + ...
        par_A2*par_covE*par_A1' +...
        par_A3*par_covE*par_A2' + ...
        par_A4*par_covE*par_A3';
    covs_true.ONB(3,:,:) = ...
        par_A2*par_covE +...
        par_A3*par_covE*par_A1' + ...
        par_A4*par_covE*par_A2' +...
        par_A1*par_covE*par_A3' + ...
        par_A2*par_covE*par_A4' +...
        par_A3*par_covE*par_A1' + ...
        par_A4*par_covE*par_A2';
    covs_true.ONB(4,:,:) = ...
        par_A3*par_covE + ...
        par_A4*par_covE*par_A1' +...
        par_A1*par_covE*par_A2' + ...
        par_A2*par_covE*par_A3' +...
        par_A3*par_covE*par_A4' + ...
        par_A4*par_covE*par_A1';
    covs_true.ONB(5,:,:) = ...
        par_A4*par_covE +...
        par_A1*par_covE*par_A1' + ...
        par_A2*par_covE*par_A2' +...
        par_A3*par_covE*par_A3' + ...
        par_A4*par_covE*par_A4';
    covs_true.ONB(6,:,:) = ...
        par_A1*par_covE + ...
        par_A2*par_covE*par_A1' +...
        par_A3*par_covE*par_A2' + ...
        par_A4*par_covE*par_A3';
    covs_true.ONB(7,:,:) = ...
        par_A2*par_covE +...
        par_A3*par_covE*par_A1' + ...
        par_A4*par_covE*par_A2';
    covs_true.ONB(8,:,:) = ...
        par_A3*par_covE + ...
        par_A4*par_covE*par_A1';
    covs_true.ONB(9,:,:) = ...
        par_A4*par_covE;
    
end

% calculate the norms
for ii = 1:(covs_true.numOfLags + 1 )   
    covs_true.ONB_norms(ii) = norm( squeeze( covs_true.ONB(ii,:,:) ));
end

% specifically save lag 0 - shortcut
covs_true.lag0 = squeeze(covs_true.ONB(1,:,:));

% recommended cut-of-lag
covs_true.cutOfLag = MA_order + 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%        set up 'sigma'

signal = trace(squeeze(covs_true.ONB(1,:,:)));
sigma = sqrt(signal/20);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              true spectral density
% pre-alocate space and  prepare frequency grid
specDensity_true = [];
specDensity_true.nGridFreq = nGridFreq;
specDensity_true.gridFreq = linspace(-pi,pi,specDensity_true.nGridFreq);

% save the spectral density
specDensity_true.ONB = zeros(nGridFreq, onb.nBasis, onb.nBasis); % this is in the ONB
if (MA_order == 2)
    for k=1:specDensity_true.nGridFreq
        omega = specDensity_true.gridFreq(k);
        specDensity_true.ONB(k,:,:) =  1/(2*pi) *...
            (eye(onb.nBasis) + par_A1*exp(-1i*omega) + par_A2*exp(-2*1i*omega))*...
            par_covE*...
            (eye(onb.nBasis) + par_A1*exp(-1i*omega) + par_A2*exp(-2*1i*omega))';
    end
end
if (MA_order == 4)
    for k=1:specDensity_true.nGridFreq
        omega = specDensity_true.gridFreq(k);
        specDensity_true.ONB(k,:,:) =  1/(2*pi) *...
            (eye(onb.nBasis)+...
            par_A1*exp(-1i*omega)   + par_A2*exp(-2*1i*omega)+...
            par_A3*exp(-3*1i*omega) + par_A4*exp(-4*1i*omega) )*...
            par_covE*...
            (eye(onb.nBasis)+...
            par_A1*exp(-1i*omega)   + par_A2*exp(-2*1i*omega)+...
            par_A3*exp(-3*1i*omega) + par_A4*exp(-4*1i*omega) )';
    end
end
if (MA_order == 8)
    for k=1:specDensity_true.nGridFreq
        omega = specDensity_true.gridFreq(k);
        specDensity_true.ONB(k,:,:) =  1/(2*pi) *...
            (eye(onb.nBasis)+...
            par_A1*exp(-1i*omega)   + par_A2*exp(-2*1i*omega)+...
            par_A3*exp(-3*1i*omega) + par_A4*exp(-4*1i*omega) +...
            par_A1*exp(-5*1i*omega) + par_A2*exp(-6*1i*omega) +...
            par_A3*exp(-7*1i*omega) + par_A4*exp(-8*1i*omega))*...
            par_covE*...
            (eye(onb.nBasis)+...
            par_A1*exp(-1i*omega)   + par_A2*exp(-2*1i*omega)+...
            par_A3*exp(-3*1i*omega) + par_A4*exp(-4*1i*omega) +...
            par_A1*exp(-5*1i*omega) + par_A2*exp(-6*1i*omega) +...
            par_A3*exp(-7*1i*omega) + par_A4*exp(-8*1i*omega))';
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              simulate the process  -  as functional process

% pre-alocate the censor struct
censor = [];

% simulate the noise process
nGridTime_noise = nGridTime + 2*dataFull_padding + 10;
noise = mvnrnd( zeros(nGridTime_noise,onb.nBasis) , par_covE );

% get the MA process
dataFullONB_padded = zeros( nGridTime + 2*dataFull_padding, onb.nBasis);
if (MA_order == 2)
    for t = 1:(nGridTime + 2*dataFull_padding)
        dataFullONB_padded(t,:) = mu_true + noise(10+t,:)' +...
            par_A1*noise(10+t-1,:)' + par_A2*noise(10+t-2,:)';
    end
end
if (MA_order == 4)
    for t = 1:(nGridTime + 2*dataFull_padding)
        dataFullONB_padded(t,:) = mu_true + noise(10+t,:)' +...
            par_A1*noise(10+t-1,:)' + par_A2*noise(10+t-2,:)' +...
            par_A3*noise(10+t-3,:)' + par_A4*noise(10+t-4,:)';
    end
end
if (MA_order == 8)
    for t = 1:(nGridTime + 2*dataFull_padding)
        dataFullONB_padded(t,:) = mu_true + noise(10+t,:)' +...
            par_A1*noise(10+t-1,:)' + par_A2*noise(10+t-2,:)' +...
            par_A3*noise(10+t-3,:)' + par_A4*noise(10+t-4,:)'+...
            par_A1*noise(10+t-5,:)' + par_A2*noise(10+t-6,:)'+...
            par_A3*noise(10+t-7,:)' + par_A4*noise(10+t-8,:)';
    end
end

% save the full data
dataFull_padded = dataFullONB_padded*onb.onbMatrix';
dataFull = dataFull_padded( (dataFull_padding+1):(nGridTime+dataFull_padding) , :) ;
dataFullONB = dataFullONB_padded(dataFull_padding + (1:nGridTime),:);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                 censored version

% create the "censor" structure - where all the sparsified data is saved
censor.nGridTime = nGridTime;
censor.onb.nGridSpace = onb.nGridSpace;
censor.onb.gridSpace = onb.gridSpace;

% censorship
censor.zeroone = zeros(onb.nGridSpace, nGridTime); % censorBool(x,n) = 1 iff the I do observe position "x" on the "n"-th curve
censor.Hspace = cell(nGridTime,1);
censor.Honb = cell(nGridTime,1);
censor.grid = cell(nGridTime,1);
censor.nGrid = zeros(nGridTime,1);
censor.dataWOnoise = cell(nGridTime,1);
censor.data = cell(nGridTime,1);
for n=1:nGridTime

    % random number of points
    nChoose = nRandi_min+randi(nRandi_max-nRandi_min); % uniform distribution 1-12
    censor.zeroone(:,n) = logical( randperm(onb.nGridSpace) <= nChoose );

    % create actual censorH matrix
    numObs = sum(censor.zeroone(:,n));
    censor.nGrid(n,1) = sum(censor.zeroone(:,n));
    censor.Hspace_act = zeros(numObs, onb.nGridSpace );
    censor.grid_act = zeros(numObs,1);
    j = 1;
    for x=1:onb.nGridSpace        
        if (censor.zeroone(x,n) == 1)
            censor.Hspace_act(j,x) = 1;
            censor.grid_act(j) = x;
            j = j+1;
        end        
    end
    censor.Hspace{n,:} = censor.Hspace_act;
    censor.Honb{n,:} = censor.Hspace_act*onb.onbMatrix;
    censor.grid{n,:} = censor.grid_act;
    censor.dataWOnoise{n,:} = censor.Hspace{n,:}*dataFull(n,:)';
    censor.data{n,:} = censor.dataWOnoise{n,:} + sigma*randn(numObs,1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       split the data, prepare for the K-FOLD cross-validation
censor.cv_K = cv_K;

% NEW - assign the same CV batch to entire curve
censor.cv_batch = cell(nGridTime,1);
for t=1:nGridTime
    censor.cv_batch{t} = ones( 1, censor.nGrid(t)) * randi(censor.cv_K);
end

% calculate how much data points are in individual batches
censor.cv_batch_counts = zeros(1,censor.cv_K);
for kk = 1:censor.cv_K
    for t=1:censor.nGridTime
        censor.cv_batch_counts(kk) = censor.cv_batch_counts(kk) + sum(censor.cv_batch{t} == kk);
    end
end




end