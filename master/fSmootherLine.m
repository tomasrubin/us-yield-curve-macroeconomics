function mu_est = fSmootherLine(X,Y,bw,degree,gridSmoother)
% Performs the line smoother of pairs (X,Y)
%
% X ... the positions of the points on the [0,1] grid
% Y ... the value
% bw ... the bandwidth parameter to be used
% degree ... what smoother to use, 1= local-linear, 2= local-quadratic
% gridSmoother ... grid of [0,1] where the smoother is evaluated


nGridSmoother = length(gridSmoother);

% smoother for mu
mu_est = zeros(nGridSmoother,1);
for ii = 1:nGridSmoother
    indx = find(abs( X(:) - gridSmoother(ii) ) <= bw + 10^(-6)); % find those that are distant at most bw (plus something small)
    deviation = X(indx,1) - gridSmoother(ii); % deviation from the grid_smoother, i.e. distance with a sign +/-
    W = spdiags( 3/4*( 1-(deviation/bw).^2 ) , 0, length(indx),length(indx)); % spatial weights
    
    % create the design matrix
    desM = zeros(length(indx),degree+1);
    desM(:,1) = 1;
    for jj = 1:degree
        desM(:,jj+1) = deviation.^jj;
    end
    
    % calculate the smoother
    beta=pinv(desM'*W*desM) * (desM'*W* Y(indx) );
    mu_est(ii) = beta(1);
end

end