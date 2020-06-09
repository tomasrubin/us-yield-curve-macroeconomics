function [censor, dataFull, dataFullONB, dataFull_padded, dataFullONB_padded, covs_true, mu_true, specDensity_true,maxEigenvalue,sigma] = fSimulate_AR2020_v3(onb,simulated_process,nGridTime,nRandi_max,nRandi_min,nGridFreq,dataFull_padding)

% some constants
nBurnIn = 10000; % burn-in for simulation of AR, in order to start from the stationary distribution
cv_K = 10; % the number of partitions for K-fold CV

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              define the sd 'sigma' 


switch simulated_process
    case 0 % IID
        maxEigenvalue = 0;
    
    case 1 % FAR 0.7        
        maxEigenvalue = 0.7;
        
    case 2 % FAR 0.9        
        maxEigenvalue = 0.9;
           
    
    otherwise % unknown process => stop the script
        fprintf('Error: unknown process\n')
        return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              define the parameters of the AR process
% autoregressive operator

g = onb.gridSpace; % short-cut for the spatial grid
parA0_inSpace =  sin((g'-g) );
parA0_inSpace = parA0_inSpace .* maxEigenvalue / max(abs(eig(parA0_inSpace)))*onb.nGridSpace; % make sure it has the designated maxEigenvalue
parA0 = onb.onbMatrix'*parA0_inSpace*onb.onbMatrix / onb.nGridSpace^2; % express parA0 in the onb.basis

% the covariance matrix of the stochastic error
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
parS0 = onb.onbMatrix'*stochCov_inSpace*onb.onbMatrix / onb.nGridSpace^2;

%% zero mean function
parMU_inSpace = zeros(1,onb.nGridSpace); % in the space variable
parMU = onb.onbMatrix' * parMU_inSpace' / onb.nGridSpace; % in the ONB
mu_true = parMU; % output in ONB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              true auto-covariances

% pre-alocate space
covs_true = [];

% calculate 
covs_true.lag0 = zeros(onb.nBasis);
for ii = 0:1000 % formula by Bosq, it should be infinite sum
    covs_true.lag0 = covs_true.lag0 + parA0^ii * parS0 * (parA0'^ii);
end

% save how many lags I have and pre-alocate space
covs_true.numOfLags = nGridTime - 1;
covs_true.ONB = zeros(covs_true.numOfLags+1, onb.nBasis, onb.nBasis);
covs_true.ONB_norms = zeros(covs_true.numOfLags+1,1);

% ssave (again) the first covariance
covs_true.ONB(1,:,:) = covs_true.lag0;
covs_true.ONB_norms(1) = norm(covs_true.lag0);

% other covariances
for ii = 2:(covs_true.numOfLags + 1 )
    covs_true.ONB(ii,:,:) = parA0* squeeze(covs_true.ONB(ii-1,:,:));        
    covs_true.ONB_norms(ii) = norm( squeeze( covs_true.ONB(ii,:,:) ));
end

% recommended cut-of-lag
covs_true.cutOfLag = find( covs_true.ONB_norms < 0.001 , 1);
if (~isnumeric(covs_true.cutOfLag) || isempty(covs_true.cutOfLag)), covs_true.cutOfLag = 100; end
if (covs_true.cutOfLag > 101), covs_true.cutOfLag = 100; end
if (covs_true.cutOfLag < 10), covs_true.cutOfLag = 10; end


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

% save the TRUE spec density here
specDensity_true.density = zeros(specDensity_true.nGridFreq, onb.nBasis, onb.nBasis); % this is in the ONB
for k = 1:specDensity_true.nGridFreq
    omega = specDensity_true.gridFreq(k);    
    specDensity_true.ONB(k,:,:) =  1/(2*pi)* ((eye(onb.nBasis)-parA0*exp(-1i*omega)) \ ... % inv(A)*b = A\b
        parS0*inv(eye(onb.nBasis)-parA0*exp(-1i*omega))');     
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%              simulate the process  -  as functional process

censor = [];


% pre-save the noise procceess
noise = mvnrnd(zeros(nBurnIn,onb.nBasis), parS0);

% burn-in, to achieve the stationary initial distribution
trajectoryONB = zeros(onb.nBasis,1);
for t = 1:nBurnIn    
    % deterministic term
    trajectoryONB = parA0 * trajectoryONB;    
    % stochastic error
    trajectoryONB = trajectoryONB + noise(t,:)';
end

% simulate the trajectory - without the common mean
nGridTime_padded = nGridTime + 2*dataFull_padding;
noise = mvnrnd(zeros(nGridTime_padded,onb.nBasis), parS0);
dataFullONB_padded = zeros( nGridTime_padded, onb.nBasis);
dataFullONB_padded(1,:) = trajectoryONB';
for t = 2:nGridTime_padded
    % deterministic term
    dataFullONB_padded(t,:) = parA0 * dataFullONB_padded(t-1,:)';    
    % stochastic error
    dataFullONB_padded(t,:) = dataFullONB_padded(t,:)' + noise(t,:)';
end

% add the common mean function
for t = 1:nGridTime_padded
    dataFullONB_padded(t,:) = dataFullONB_padded(t,:) + parMU';
end

% express in spatial variable
dataFull_padded = dataFullONB_padded*onb.onbMatrix';
dataFull = dataFull_padded(dataFull_padding + (1:nGridTime),:);
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