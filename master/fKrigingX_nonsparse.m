function [kriging, kriging_padded] = fKrigingX_nonsparse(censor,onb,covLags,mu_est_ONB,sigma_est,dynamic,calculate_var, dataFull_padding)
% This function performs the "functional data recovery" a.k.a. kriging
%
% censor ... the data of the process X
% onb ... the ONB struct
% covLags ... struct of either true or estimated auto-covariances
% calculate_var ... 1 = returns also the variance for each functional datum
%                   0 = the var fields will be set NAN
% dataFull_padding ... how much to pad at both ends of the sample

% how much time is in the padded sample
nGridTime_padded = censor.nGridTime + 2*dataFull_padding;

% pre-alocate the space
kriging_padded = [];
kriging_padded.est = nan(nGridTime_padded, onb.nBasis);
kriging_padded.var = nan(nGridTime_padded, onb.nBasis, onb.nBasis);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dynamic %%          dynamic kriging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    %% prepare the prior matrix for the functional data (in onb.basis)
    dimX = onb.nBasis * nGridTime_padded; % vectorize the functional time series
    prior_covX = zeros(dimX);        
    for lag = 0:covLags.cutOfLag
        for t = 1:(nGridTime_padded-lag)
            s = t+lag;
            % prepare the covariance between X_t and X_s        
            covAct = squeeze(covLags.ONB(lag+1,:,:)); 
            for ii = 1:onb.nBasis            
                for jj = 1:onb.nBasis
                    % fill in the covariance between < X_t, e_ii > and <X_s, e_jj>
                    % i.e. elements: (t-1)*onb.nBasis + ii .... (s-1)*onb.nBasis + jj
                    prior_covX( (t-1)*onb.nBasis + ii, (s-1)*onb.nBasis +jj) = covAct(jj,ii); %covAct(ii,jj);
                    prior_covX( (s-1)*onb.nBasis + ii, (t-1)*onb.nBasis +jj) = covAct(ii,jj); %covAct(jj,ii);
                end
            end
        end
    end

    % save it as sparse matrix
    prior_covX = sparse(prior_covX);
    % nnz(prior_covX) / prod(size(prior_covX)) % the achieved sparsity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% prepare the prior for mu
    prior_muX = zeros(1,dimX);
    for t = 1:nGridTime_padded
        for ii = 1:onb.nBasis
            prior_muX( (t-1)*onb.nBasis + ii ) = mu_est_ONB(ii);
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% prepare the mapping operator to the functional data (in onb.basis)
    dimZ = sum(censor.nGrid); % the number of the observables
    obsVec = zeros(dimZ,1); % vectorized data - first few elements correspond to the first censored X_1
    start_skip_Z = zeros(censor.nGridTime,1); % the element on the position "t" says that before saving censor.data{t} I
    for t=1:censor.nGridTime 
        if t>1, start_skip_Z(t) = sum( censor.nGrid(1:(t-1)) ); end
        for ii = 1:censor.nGrid(t)        
            obsVec( start_skip_Z(t) + ii ) = censor.data{t}(ii);
        end
    end

    obsOperator = zeros(dimZ, dimX);
    % obsOperator(ii,jj) - how the jj-th onb.basis member of X contributes to the ii-th observable
    for t = 1:censor.nGridTime
        for ii = 1:onb.nBasis
            H_act = censor.Honb{t};
            for jj = 1:censor.nGrid(t)
                obsOperator( start_skip_Z(t)+jj, dataFull_padding*onb.nBasis + (t-1)*onb.nBasis + ii ) = H_act(jj,ii);
            end        
        end
    end


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%     Gaussian completion

    % this is the model: obsVec = obsOperator * X + error
    B = sparse(obsOperator);

    % kriging
    posterior_muX = prior_muX' + prior_covX * B'*((B*prior_covX*B'+sigma_est^2*eye(dimZ)) \ (obsVec - B*prior_muX'  )); % INV(A)*b == A\b
    if calculate_var
        posterior_covX = prior_covX-prior_covX*B'* ((B*prior_covX*B'+ sigma_est^2*eye(dimZ)) \ (B*prior_covX));        
    end

    % unvectorize posterior mean
    for t=1:nGridTime_padded
        for ii=1:onb.nBasis
            kriging_padded.est( t,ii ) = posterior_muX( (t-1)*onb.nBasis + ii );        
        end
    end

    % unvectorize the variance
    if calculate_var
        for t=1:nGridTime_padded
            for ii=1:onb.nBasis
                for jj=1:onb.nBasis
                    kriging_padded.var(t,ii,jj) = posterior_covX(  (t-1)*onb.nBasis + ii,  (t-1)*onb.nBasis + jj );
                end
            end
        end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
else %%       static kriging
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % prepare prior covariance    
    prior_cov = covLags.lag0;
    [V,D] = eig(prior_cov);
    diagD = diag(D);
    diagD = max(diagD,1e-6);
    prior_cov = V * diag(diagD) * inv(V); % make it positive definite
    
        
    % cycle for kriging
    for t = 1:nGridTime_padded
        
        % short-cut for the censor operator
        A = censor.Honb{t};
        
        % kriging of the mean        
        kriging_padded.est(t,:) = mu_est_ONB + prior_cov*A'*inv(A*prior_cov*A'+sigma_est^2*eye(censor.nGrid(t))) * (censor.data{t} - censor.Honb{t}*mu_est_ONB);
        
        % kriging of the var
        if calculate_var
            kriging_padded.var(t,:,:) = prior_cov - prior_cov*A'*inv(A*prior_cov*A'+sigma_est^2*eye(censor.nGrid(t)))*A*prior_cov;
        end
        
    end    

end


%% non-padded version
kriging = [];
kriging.est = kriging_padded.est( (dataFull_padding+1):(dataFull_padding+censor.nGridTime) ,:);
kriging.var = kriging_padded.var( (dataFull_padding+1):(dataFull_padding+censor.nGridTime) ,:,:);

end