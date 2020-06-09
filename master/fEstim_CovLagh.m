function covLagh_smoother = fEstim_CovLagh(censor,onb,lag,mu_est_inSpace,bw_r)

% performs a smoother for lag_h (!= 0) auto-covariance
if lag == 0
    fprintf('Error: lag == 0\n')
    return;
end

% if the lag_h is negative 
transpose_at_the_end = 0;
if lag < 0
    lag = -lag;
    transpose_at_the_end = 1;
end

% if lag = 0, remove diagonal
remove_diag = 0;
if lag == 0, remove_diag = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% here go all raw covariances


% calculate the number of lag-h covariances

% calculate the number of raw covariances
nRawCovs = 0;
for t = 1:(censor.nGridTime-lag)
    if remove_diag % lag == 0
        nRawCovs = nRawCovs + censor.nGrid(t)^2-censor.nGrid(t);    
    else % lag > 0
        nRawCovs = nRawCovs + censor.nGrid(t+lag)*censor.nGrid(t);    
    end    
end

% get raw covariances (without the diagonal)
iRawCov = 1; % the index of the raw cow now being filled
X = zeros(nRawCovs,2); % the "covariates"
Y = zeros(nRawCovs,1); % the value of the raw covariates product
for t = 1:(censor.nGridTime-lag)
    for ii = 1:censor.nGrid(t+lag)
        for jj = 1:censor.nGrid(t)
            if (remove_diag == 0) || ( ii ~= jj) % include it iff it's not the diag of lag 0
                X(iRawCov,1) = onb.gridSpace(censor.grid{t+lag}(ii));
                X(iRawCov,2) = onb.gridSpace(censor.grid{t}(jj));
                Y(iRawCov) = (censor.data{t+lag}(ii) - mu_est_inSpace(censor.grid{t+lag}(ii)) )*...
                    (censor.data{t}(jj) - mu_est_inSpace(censor.grid{t}(jj)) );
                iRawCov = iRawCov+1;    
            end
        end
    end
end

%
covLagh_smoother = zeros(onb.nGridSmoother,onb.nGridSmoother);
for ii = 1:onb.nGridSmoother
    for jj = 1:onb.nGridSmoother
        % now we're filling up the point covLag0(ii,jj), i.e. [onb.gridSpace(ii), onb.gridSpace(jj)]                
        % choose only those indeces where I'm less then "bw" far away
        indx = find( abs(X(:,1) - onb.gridSmoother(ii)) <= bw_r + 10^(-6) & abs( X(:,2) - onb.gridSmoother(jj)) <= bw_r + 10^(-6) );
        dist1 = X(indx,1) - onb.gridSmoother(ii);        
        dist2 = X(indx,2) - onb.gridSmoother(jj);        
        % define the weights 
        W = spdiags( 9/16*(1-(dist1/bw_r).^2).*(1-(dist2/bw_r).^2), 0, length(indx),length(indx) );
%         W = diag( 9/16*(1-(dist1/bw_r).^2).*(1-(dist2/bw_r).^2) );
        % get the design matrix
        desM = zeros( length(indx), 3 );
        desM(:,1) = 1;
        desM(:,2) = dist1;
        desM(:,3) = dist2;
        beta=pinv(desM'*W*desM) * (desM'*W* Y(indx) );
        covLagh_smoother(ii,jj)=beta(1);
    end
end


% if it was a negative lag, transpose
if transpose_at_the_end
    covLagh_smoother = covLagh_smoother';
end



end