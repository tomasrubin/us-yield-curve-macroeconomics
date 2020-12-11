%% A case study analysing the dependence of the US Treasury yield curve on the macroeconomic variables.
% 
% This script runs the analysis of the US Treasury yield curve dependence on the macroeconomic variables
% in the US economy: the annual change in industrial production, the annual inflation, and the federal
% funds rate. The methodology follows the novel non-parametric toolbox for spectral analysis of sparsely
% observed functional time series. For the details of the methodology and the analysis, see [1]. The code
% for this case study is openly available on GitHub [2].
%
% Individuals are free to use the codes for the purpose academic research, provided it is properly acknowledged.
% For any other use, permission must first be arranged with the author. Unless otherwise specified, the author
% Please contact me if you find errors in the codes.
%
% Running the script: press F5 to run the entire file. The results are outputed as plots.
%
% Author: Tomas Rubin
%         tomas.rubin@gmail.com
%         linkedin.com/in/tomas-rubin/
%         github.com/tomasrubin/
%
% References:
%   [1] Rubin, T. (2020). "US Treasury yield curve and macroeconomy interaction: evidence from the non-parametric functional lagged regression approach." arXiv preprint
%   [2] github.com/tomasrubin/us-yield-curve-macroeconomics
% 
% 
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  preliminaries

% here's the library
addpath('master')
addpath('master/fdaM')

% prepare the ONB basis
onb = fPrepareONB(30*12+1, 12, 30); % each year has 12 grid points;

% the resolution of the frequency grid on [-pi,pi]
nGridFreq = 1000;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  load data

macro_table_full = readtable('us_macro.xlsx');
data_table_full = readtable('us_yields.xlsx');

data_table_full.year = year(data_table_full.Date);
data_table_full.month = month(data_table_full.Date);
data_table_full.day = day(data_table_full.Date);


%% extract macro data
begin_year = 1985;
end_year = 2000;

macro = [];
macro.date = macro_table_full.Date( (year(macro_table_full.Date) >= begin_year) & (year(macro_table_full.Date) <= end_year) );
macro.annual_industry = macro_table_full.annualIndustryGrowth( (year(macro_table_full.Date) >= begin_year) & (year(macro_table_full.Date) <= end_year) );
macro.inflation = macro_table_full.annualInflation( (year(macro_table_full.Date) >= begin_year) & (year(macro_table_full.Date) <= end_year) );
macro.rate = macro_table_full.targer_rate( (year(macro_table_full.Date) >= begin_year) & (year(macro_table_full.Date) <= end_year) );
macro.x = [100*macro.annual_industry, 100*macro.inflation, macro.rate];
macro.nGridTime = length(macro.date);
macro.year  = year( macro.date );
macro.month = month( macro.date );
macro.year_month = year( macro.date  ) * 100 + month( macro.date  );



%% extract the yields corresponding to the months of macro variables
% find yield in the middle of the month

data = [];
data.year_month = macro.year_month;
data.nGridTime = length(data.year_month);

% maturities
data.maturities_gs_indx = 1+[1,2,3, 6, 12, 2*12, 3*12, 5*12, 7*12, 10*12, 20*12, 30*12]; % indecies on the spatial grid
data.maturities_real = [1/12,2/12,3/12, 6/12, 1, 2, 3, 5, 7, 10, 20, 30]; % point on [0,30]

% yields
data.y = nan(data.nGridTime,12);
for t = 1:data.nGridTime
    
    % extract the current month, from day 15
    yields_sub = data_table_full( (data_table_full.year == macro.year(t)) & (data_table_full.month == macro.month(t) ) & (data_table_full.day >= 15 ) , :);
    
    % extract the day 15, or later
    data.yields(t,:) = yields_sub{1, 2:13};
end

% keep only those maturities, for which there is at least one observation
keep = (sum(~isnan(data.yields),1)>0);
keep(11) = 0;
data.yields = data.yields(:,keep);

% differences
% data.yields(2:end,:) = data.yields(2:end,:) - data.yields(1:(data.nGridTime-1),:);
% data.yields(1,:) = data.yields(2,:); % XXX this should be taken care of later

data.maturities_real = data.maturities_real(:,keep);
data.maturities_gs_indx = data.maturities_gs_indx(:,keep);
data.n_maturities = length(data.maturities_real);

% change the maturities on the gridSpace so they have equidistant spacing
data.maturities_gs_indx = ceil(linspace(1, onb.nGridSpace, data.n_maturities ));


% visiualise the transofrmation of maturities to the equidistant grid
%subplot(1,1,1)
% plot(  onb.gridSpace(data.maturities_gs_indx),data.maturities_real, '-x' )
figure(1)
x = onb.gridSpace(data.maturities_gs_indx);
y = data.maturities_real;
xx = onb.gridSpace;
yy = spline(x,y,xx);
plot(x,y,'o',xx,yy,'MarkerFaceColor', 'b')
ylabel("real maturity in [0,30]")
xlabel("transformed maturity in [0,1]")
title("transformation of the spatial domain")

data.gridSpace_to_maturities = yy;

% ONB
onb = fPrepareONB(30*12+1, 12, 30); % each year has 12 grid points;
data.onb = onb;


%% nonclassical surface plot - visualise the evolution at maturities
figure(2)
subplot(1,2,1)

year_str = cell(length(data.year_month),1);
month_str = cell(length(data.year_month),1);
for ii = 1:length(data.year_month)
    year_str(ii) = extractBetween(num2str(data.year_month(ii)),1,4);
    month_str(ii) = extractBetween(num2str(data.year_month(ii)),5,6);
end


gridTime = 1:data.nGridTime;
% first curve
surf( data.maturities_real, 1:data.nGridTime, data.yields, 'edgecolor','none', 'facealpha', 0.5)
hold on
for m = 1:data.n_maturities
    plot3( data.maturities_real(m)*ones(1,data.nGridTime), gridTime, data.yields(:,m), 'k','LineWidth', 5 )
end
hold off
xlabel("maturity [years]")
zlabel("yield [%]")
ylabel("year/month")
y_years = ((1:8)-1) * 24 + 1;
yticks( y_years )
yticklabels(year_str(y_years) + "/01" )

% nonclassical surface plot - visualise the yield curves

subplot(1,2,2)
% first curve
surf( data.maturities_real, 1:data.nGridTime, data.yields, 'edgecolor','none', 'facealpha', 0.5)

% t=1;
% plot3(  data.maturities_real(maturities_keep),t*ones(size(maturities_keep)), data.yields(t,maturities_keep), 'k', 'LineWidth', 1 )
hold on
for t = 1:1:data.nGridTime
    plot3(  data.maturities_real, t*ones(size(data.maturities_real)), data.yields(t,:), 'k', 'LineWidth', 2 )
end
hold off
xlabel("maturity [years]")
zlabel("yield [%]")
ylabel("year/month")
yticks( y_years )
yticklabels(year_str(y_years) + "/01" )

suptitle("US treasury yield curve")

% %%
% 
% for t_i=1:3
%     t = t_i * 60;
%     subplot(3,1,t_i)
%     plot( data.maturities_real, data.yields(t,:), 'xk', 'LineWidth', 2)
%     title( year_str(t) +"/"+month_str(t) )
%     xlabel("maturity [years]")
%     ylabel("yield [%]")
%     ylim([4, 10])
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% visualise macro

figure(3)
plot( macro.date, macro.x, "LineWidth", 2 )
%hline( mean(macro.x) )
legend(["annual industry change","annual inflation","federal funds rate"])
ylabel("percantage points [%]")
title("macroeconomic variables")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% time-domain analysis of macro

macro.mean = mean(macro.x);

% covariance
macro.numOfCovs = macro.nGridTime+1;
macro.covs = nan(macro.numOfCovs,3,3);
for lag_idx = 1:macro.numOfCovs
    lag_real = lag_idx - 1;
    
    for ii = 1:3
        for jj = 1:3
            macro.covs(lag_idx,ii,jj) = 1/macro.nGridTime *  sum( (macro.x((1+lag_real):end,ii) - macro.mean(ii)) .*...
                (macro.x(1:(macro.nGridTime-lag_real),jj)  - macro.mean(jj) ) );
        end
    end
    
    
end

% correlation
corr_normalizer = diag( squeeze(macro.covs(1,:,:)) );
macro.corrs = nan(macro.numOfCovs,3,3);
for lag_idx = 1:macro.numOfCovs
    lag_real = lag_idx - 1;
    macro.corrs(lag_idx,:,:) = squeeze(macro.covs(lag_idx,:,:)) ./ sqrt( corr_normalizer'.*corr_normalizer );
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% spectral analysis of macro
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
% choose weight function
% weight_function = @(x) (abs(x)<=0.5) .* ( 1-6*x.^2 + 6*abs(x).^3 ) + (abs(x)>0.5) .* (2*(1-abs(x)).^3); % parzen weight
weight_function = @(x) 1-abs(x); % bartlett weight

% bartlett's span parameter
% macro.q_bartlett = macro.nGridTime;
macro.q_bartlett = 14; % XXX
macro.spectral_weights = weight_function( ((0:macro.q_bartlett))/macro.q_bartlett );

% spectral grid
macro.nGridFreq = 1000;
macro.gridFreq = linspace(-pi,pi,macro.nGridFreq);
macro.specDensity = nan(macro.nGridFreq, 3, 3);
macro.specNorm = nan(macro.nGridFreq,1);

% loop to estimate the spec density
for k = 1:macro.nGridFreq
    omega = macro.gridFreq(k);
    
    % start from lag=0
    macro.specDensity(k,:,:) = 1/(2*pi) * squeeze( macro.covs(1,:,:) );
    
    % other lags, both positive and negative
    for lag_real = 1:macro.q_bartlett
        lag_inx = lag_real+1;
        macro.specDensity(k,:,:) = squeeze( macro.specDensity(k,:,:) ) + ...   
            macro.spectral_weights(lag_inx) *...
            1/(2*pi) * ( squeeze( macro.covs(lag_inx,:,:) ) * exp(-1i*omega*lag_real) + squeeze( macro.covs(lag_inx,:,:) )' * exp(+1i*omega*lag_real) );
    end
    
    % norm of spec density
    macro.specNorm(k) = norm(squeeze(macro.specDensity(k,:,:)), 'fro');
end

macro.specDensity_positified = zeros(nGridFreq,3,3);
macro.specDensity_harmonic_eigenvalues = zeros(nGridFreq,3);
% positify the spectral density, from the midle to the right
for k = ceil(macro.nGridFreq/2):macro.nGridFreq
    
    % positify the spectral density
    m = squeeze(macro.specDensity(k,:,:));
    [V,D] = eig(m); % i.e. the decomposition:  m = V*D*V'
    diagD = real(diag(D)); % note that the harmonic eigenvalues must be always real
    macro.specDensity_harmonic_eigenvalues(k,:) = diagD;
    diagD( diagD< 0 ) = 0; % truncate those elements that are less than (numerical) zero
    m = V*diag(diagD)*V';
    macro.specDensity_positified(k,:,:) = (m+m')/2; % again, make sure it's self-adjoint
    
end
% positify the spectral density, from the midle to the left
for k_indx = 1:floor(macro.nGridFreq/2)
    k = floor(macro.nGridFreq/2) - k_indx + 1;
    
    % positify the spectral density
    m = squeeze(macro.specDensity(k,:,:));
    [V,D] = eig(m); % i.e. the decomposition:  m = V*D*V'
    diagD = real(diag(D)); % note that the harmonic eigenvalues must be always real
    macro.specDensity_harmonic_eigenvalues(k,:) = diagD;
    diagD( diagD< 0 ) = 0; % truncate those elements that are less than (numerical) zero
    m = V*diag(diagD)*V';
    macro.specDensity_positified(k,:,:) = (m+m')/2; % again, make sure it's self-adjoint
    
end



% visualise the spectral density
k_to_visualise = [500 520 540 560 600 800];

real_minmax = [min(real(macro.specDensity(:))), max(real(macro.specDensity(:)))];
imag_minmax = [min(imag(macro.specDensity(:))), max(imag(macro.specDensity(:)))];


for ii = 1:length(k_to_visualise)
    k = k_to_visualise(ii);
    omega = macro.gridFreq(k); % real omega: -pi to pi
    
    m = squeeze(macro.specDensity(k,:,:));
    
    subplot(3,length(k_to_visualise),ii)
    surf( real(m), 'EdgeColor','none')
    title("omega = "+num2str(round(omega,2)))
    if ii == 1, zlabel("real part"); end
    zlim(real_minmax)
    
    subplot(3,length(k_to_visualise),ii+length(k_to_visualise))
    surf( imag(m), 'EdgeColor','none')
    if ii == 1, zlabel("imag part"); end
    zlim(imag_minmax)
    
end

% draw norms over frequencies
subplot(3,1,3)
plot( macro.gridFreq, macro.specDensity(:,1,1) )
hold on
plot( macro.gridFreq, macro.specDensity(:,2,2) )
plot( macro.gridFreq, macro.specDensity(:,3,3) )
hold off


ylabel("marginal spec densities")
xlabel("frequency (highlighted the frequencies visualised above)")
vline( macro.gridFreq(k_to_visualise) )
suptitle("estimated spectral density \{F_\omega^X\}")

legend(["annual industry change","annual inflation","federal funds rate"])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% fit the mean function for YIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5)

bw_mu = 0.1;
mu_est = fSmootherMU_grid_unif(data, bw_mu);

% plot with equidistant maturities
subplot(2,1,1)
plot( repmat( onb.gridSpace(data.maturities_gs_indx),data.nGridTime,1) + 0.006*randn(size(data.yields)) , data.yields ,'.')
hold on
plot( onb.gridSpace, mu_est.inSpace, 'k', 'LineWidth',4)
hold off
title("mu estimated, equidistant spacing of maturities")
ylabel("yield [%]")
xlabel("maturity [years]")

% plot with real maturities
subplot(2,1,2)
plot( repmat( data.maturities_real,data.nGridTime,1) + 0.1*randn(size(data.yields)) , data.yields ,'.')
hold on
plot( data.gridSpace_to_maturities, mu_est.inSpace, 'k', 'LineWidth',4)
hold off
title("mu estimated, actual spacing of maturities")
ylabel("yield [%]")
xlabel("maturity [years]")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% cross spectral analysis between MACRO and YIELDS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)

% settings for the crossSpecDensity estimation
bw_f = 0.1;
nGridFreq = 1000;
q_bartlett = macro.q_bartlett; % XXX

crossSpecDensity_1 = fEstimateCrossSpecDensity_grid_unif( macro.x(:,1), data, mu_est, bw_f, nGridFreq, q_bartlett);
crossSpecDensity_2 = fEstimateCrossSpecDensity_grid_unif( macro.x(:,2), data, mu_est, bw_f, nGridFreq, q_bartlett);
crossSpecDensity_3 = fEstimateCrossSpecDensity_grid_unif( macro.x(:,3), data, mu_est, bw_f, nGridFreq, q_bartlett);

crossSpecDensity_ONB_joint = nan( macro.nGridFreq, onb.nBasis, 3 );
crossSpecDensity_ONB_joint(:,:,1) = crossSpecDensity_1.ONB;
crossSpecDensity_ONB_joint(:,:,2) = crossSpecDensity_2.ONB;
crossSpecDensity_ONB_joint(:,:,3) = crossSpecDensity_3.ONB;



% visualise the cross spectral densities

m = crossSpecDensity_1.smoother;
subplot(2,1,1)
surf( 30*onb.gridSmoother, crossSpecDensity_3.gridFreq, real(m) )
xlabel("maturity [years]")
ylabel("frequency")
zlabel("real part")
suptitle("cross-spectral density between yield curve and federal funds rate")


subplot(2,1,2)
surf( 30*onb.gridSmoother, crossSpecDensity_3.gridFreq, imag(m) )
xlabel("maturity [years]")
ylabel("frequency")
zlabel("imaginary part")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% recover the spectral transfer function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(7)

% no regularisation is required
Tikhonov_regularisation = 0.000;

macro.specDensity_inv = nan( macro.nGridFreq, 3, 3 );
specTransfer_ONB_joint = nan( macro.nGridFreq, onb.nBasis, 3 );
specTransfer_ONB_1 = nan( macro.nGridFreq, onb.nBasis );
specTransfer_ONB_2 = nan( macro.nGridFreq, onb.nBasis );
specTransfer_ONB_3 = nan( macro.nGridFreq, onb.nBasis );
for k = 1:macro.nGridFreq
        
    % TIKHONOV REGULARISATION
    m = squeeze(macro.specDensity_positified(k,:,:));
    m_inv = inv(m + Tikhonov_regularisation * eye(3));
    
    % save the inverted spectral density
    macro.specDensity_inv(k,:,:) = m_inv;
    
    % recover the cross spectral density
    specTransfer_ONB_joint(k,:,:) = squeeze(crossSpecDensity_ONB_joint(k,:,:)) * m_inv;
    
    specTransfer_ONB_1(k,:) = specTransfer_ONB_joint(k,:,1);
    specTransfer_ONB_2(k,:) = specTransfer_ONB_joint(k,:,2);
    specTransfer_ONB_3(k,:) = specTransfer_ONB_joint(k,:,3);
end

% integrate back to get filter coefficients
transfer = [];
transfer.nGridFreq = macro.nGridFreq;
transfer.gridFreq = macro.gridFreq;
transfer.numOfLags = 30;

% here are the joint filter coefficients
transfer.ONB_joint = zeros(2*transfer.numOfLags+1,onb.nBasis,3);
transfer.inSpace_joint = nan(2*transfer.numOfLags+1,onb.nGridSpace,3);
for lag_indx = 1:(2*transfer.numOfLags+1)
    lag_real = lag_indx-transfer.numOfLags-1;
    
    for k=1:transfer.nGridFreq
        omega=transfer.gridFreq(k);
        transfer.ONB_joint(lag_indx,:,:) = squeeze(transfer.ONB_joint(lag_indx,:,:)) +...
            1/ transfer.nGridFreq * squeeze(specTransfer_ONB_joint(k,:,:)) * exp(1i*lag_real*omega);
    end
    
    % get the transfer functions in gridSpace
    for macro_variable_i = 1:3
        transfer.inSpace_joint(lag_indx,:,macro_variable_i) = real(onb.onbMatrix * squeeze(transfer.ONB_joint(lag_indx,:,macro_variable_i))');
    end
end



% visulise the filters
lags_to_show = fliplr( -4:4 );
macro_variables_names = {"industry","inflation","policy rate"};
for ii = 1:length(lags_to_show)
    lag_real = lags_to_show(ii);
    lag_indx = lag_real + transfer.numOfLags + 1;
    for macro_variable_i = 1:3
        
        
        % second line
        v = transfer.inSpace_joint(lag_indx,:,macro_variable_i);
        subplot(3,length(lags_to_show), ii + length(lags_to_show)*(macro_variable_i-1)*1 )
        plot( data.gridSpace_to_maturities, v ,'r','linewidth',3)
%         ylim([-1 1])
        ylim([-0.5 0.5] * ceil(2*max(max(abs(transfer.inSpace_joint(:,:,macro_variable_i))))) )
        title("lag = "+num2str(lag_real))

        % label on y axis
        if ii == 1
            ylabel({ macro_variables_names{macro_variable_i},"true maturities"})
        end
        
    end
    
end

suptitle("estimated filter coefficients for the lagged regression model")


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% prediction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(8)

predictions = [];
predictions.ONB = nan(data.nGridTime, onb.nBasis);
predictions.inSpace = nan(data.nGridTime, onb.nGridSpace);
predictions.yields = nan(data.nGridTime, data.n_maturities);

for t = 1:data.nGridTime
    
    % firstly start with the ONB
    predictions.ONB(t,:) = mu_est.ONB; % start with the mean    
    for lag_indx = 1:(2*transfer.numOfLags+1)
        lag_real = lag_indx-transfer.numOfLags-1;
        
        t_x = (t-lag_real); % this is the time of MACRO that I'm going to add now
        if (t_x >= 1) && (t_x <= data.nGridTime) % only if I have MACRO data for this "t_x"
            for macro_id = 1:3                
                predictions.ONB(t,:) = predictions.ONB(t,:) + real(transfer.ONB_joint(lag_indx,:,macro_id) * (macro.x(t_x,macro_id) - macro.mean(macro_id)));
            end            
        end        
    end    
    
    % predictions inSpace and yields
    predictions.inSpace(t,:) = onb.onbMatrix * predictions.ONB(t,:)';
    predictions.yields(t,:) = predictions.inSpace(t,data.maturities_gs_indx);
    
end


% visualise
t_visualise = sort(randsample(data.nGridTime,12));

vis_lim = [ floor(min(min(predictions.yields(:)),min(data.yields(:)))) ,  ceil(max(max(predictions.yields(:)),max(data.yields(:)))) ];

subplot(1,1,1)
for ii = 1:length(t_visualise)
    t = t_visualise(ii);
    
    subplot(3,4,ii)
    plot( data.gridSpace_to_maturities, predictions.inSpace(t,:), '-b' )
    hold on
    scatter( data.maturities_real, predictions.yields(t,:), 'filled', 'ob' )
    scatter( data.maturities_real, data.yields(t,:), 'filled', 'or' )
    hold off
    ylim(vis_lim)
    title( year_str(t) +"/"+month_str(t) )
    
    
end

suptitle("predicted (blue) and true (red) yields")



%% asses the prediction performance
m = repmat( mean(data.yields), data.nGridTime, 1 );

mse_total = mean(( m(:) - data.yields(:)).^2 );
mse_residual = mean(( predictions.yields(:) - data.yields(:)).^2 );


r2 = 1 - mse_residual / mse_total 



















