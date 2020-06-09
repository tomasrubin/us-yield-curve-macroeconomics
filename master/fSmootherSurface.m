function cov_est = fSmootherSurface(X,Y,bw,gridSmoother,symmetric)
% Performs the local-linear surface smoother on tripplets (X(:,1),X(:,2),Y)
%
% X ... the positions of the points on the [0,1] grid
% ....... X(:,1) the first coordinate
% ....... X(:,2) the second coordinate
% Y ... the value
% bw ... the bandwidth parameter to be used
% gridSmoother ... grid of [0,1] where the smoother is evaluated
% symmetric ... logical, if the surface to be estimated is symmetric, then
% ... only one triangle is to be evaluated

nGridSmoother = length(gridSmoother);


% smoother for cov_est
cov_est = zeros(nGridSmoother,nGridSmoother);
for ii = 1:nGridSmoother
    
    % take into account anticipated symmetry
    jj_start = 1;
    if symmetric, jj_start = ii; end
    
    for jj = jj_start:nGridSmoother
        % now we're filling up the point covLag0(ii,jj), i.e. [gridSpace(ii), gridSpace(jj)]                
        % choose only those indeces where I'm less then "bw" far away
        indx = find( abs(X(:,1) - gridSmoother(ii)) <= bw + 10^(-6) & abs( X(:,2) - gridSmoother(jj)) <= bw + 10^(-6) );
        dist1 = X(indx,1) - gridSmoother(ii);        
        dist2 = X(indx,2) - gridSmoother(jj);        
        % define the weights 
        W = spdiags( 9/16*(1-(dist1/bw).^2).*(1-(dist2/bw).^2), 0, length(indx),length(indx) );
        % get the design matrix
        desM = zeros( length(indx), 3 );
        desM(:,1) = 1;
        desM(:,2) = dist1;
        desM(:,3) = dist2;
        beta=pinv(desM'*W*desM) * (desM'*W* Y(indx) );
        cov_est(ii,jj)=beta(1);
    end
end

% complete the matrix for symetric case
if symmetric
    for ii = 2:nGridSmoother
        for jj = 1:(ii-1)
            cov_est(ii,jj) = cov_est(jj,ii);
        end
    end
end

end