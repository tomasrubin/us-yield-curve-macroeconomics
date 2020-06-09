function censor = fCreate_censor( data )


censor = [];

% shortcut for ONB
onb = data.onb;

% create the "censor" structure - where all the sparsified data is saved
censor.nGridTime = data.nGridTime;
censor.onb.nGridSpace = onb.nGridSpace;
censor.onb.gridSpace = onb.gridSpace;

% censorship
% censor.zeroone = zeros(onb.nGridSpace, nGridTime); % censorBool(x,n) = 1 iff the I do observe position "x" on the "n"-th curve
censor.Hspace = cell(data.nGridTime,1);
censor.Honb = cell(data.nGridTime,1);
censor.grid = cell(data.nGridTime,1);
censor.nGrid = zeros(data.nGridTime,1);
censor.data = cell(data.nGridTime,1);

for t=1:data.nGridTime

    % identify what maturities are present at time "t"
    maturities_present = ~isnan(data.yields(t,:));
    
    % fill in the censored data
    numObs = sum(maturities_present);
    censor.nGrid(t) = numObs;
    censor.grid{t} = data.maturities_gs_indx(maturities_present);
    censor.data{t} = data.yields(t, maturities_present);
    
    % the censor matrices
    Hspace_act = zeros(numObs, onb.nGridSpace);
    j = 1;
    for x=1:onb.nGridSpace
        if any( x== censor.grid{t} )
            Hspace_act(j,x) = 1;
            j = j+1;
        end
    end
    censor.Hspace{t} = Hspace_act;
    censor.Honb{t} = Hspace_act*onb.onbMatrix;
    
end