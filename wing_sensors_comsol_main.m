% This script finds optimal sensor locations using the SSPOC method. Run
% this script with optimal_sensing_FEAwing as the working directory.
%
% Steps:
%   1) Load simulated strain data from COMSOL 
%   2) Transform simulated data with neural sensors
%   3) Dimensionality reduction on sensor data (PCA)
%   4) Optimization to find sensor locations
%
% Makes use of the cvx package, available from: http://cvxr.com/cvx/

%% specify COMSOL simulation parameters and save preferences

OM      = '2OM';   % order of magnitude of stiffness gradient ('0OM' or '2OM')
axName  = 'Yaw';   % axis of rotation
rot     = '1p0';   % rotation rate, in rad/s
stiffs  = 3e9;     % average stiffness, in Pa (column vector if multiple values)
damps   = 1;       % damping ratio (column vector if multiple values)
tStep   = '2e-4';  % time step of COMSOL simultion, in s

simIter = 4;       % number of iterations to run optimization 
suff    = ['_' OM '_test'];  % suffix for filename to save results to
saveResults = 1;   % 1 to save results to file, 0 does not save
loadResults = 0;   % 1 to check for existing results file with same name, 0 to create new file or overwrite existing file
plotFlag    = 1;   % 1 to plot results of each iteration, 0 to suppress plots

homedir     = pwd;  % path to optimal_sensing_FEAwing
                    % COMSOL simulation data must be inside /simulated_wing_data/ within this directory
resultsDir  = [pwd '/results'];  % path to save results to


%% parameters for sensor optimization

Pars.sampFreq      = 1e4; % sampling frequency for sensor simulations (need not be the same as COMSOL sampling frequency)
Pars.chordElements = 26;  % chord elements on wing
Pars.spanElements  = 51;  % span elements on wing
Pars.flapFrequency = 25;  % wingbeat frequency (Hz)

Pars.simStartup = .04;  % simulation startup time to trim off, in s
Pars.simEnd     = .2;   % simulation end time, in s

% neural simulation parameters (see neuralTransormationOfData.m)
Pars.staFreq        = 1;    
Pars.staWidth       = 4; 
Pars.staDelay       = 5;
Pars.nldShift       = 1e-4; 
Pars.nldGrad        = 5e5;  
Pars.normalizeVal   = 1;    
Pars.subSamp        = 1;

Pars.spikeReps     = 100;  % how many repetitions of spikes to simulate
Pars.refPer        = 15;   % refractory period, ms 

Pars.E = stiffs;          % wing stiffness
Pars.alpha = damps*100;   % damping ratio x100

% parameters for splitting into training and test sets
Pars.trainFraction  = 0.9; % fraction of data to use for training vs test

% parameters for SSPOC
Pars.rmodes             = 20;   % how many modes to keep during dim reduction
Pars.wTrunc             = 3;    % desired number of sensors; choose an integer between 1 and 30
Pars.elasticNet         = 0.9;
Pars.penaltyVec         = 0;    % must be row vector [1, Pars.chordElements x Pars.spanElements] or 0
Pars.eps                = 1e-6; % tolerance for w reconstruction (only used if numClasses>2)
Pars.columnCoupling     = 0;    % column coupling weight to promote nonzero rows of s (only used if numClasses>2)

ParCombo = parCombo(Pars); % set up for parameter sweep
parComboTable = struct2table(ParCombo,'AsArray',1);


%% set up to save results
resultsFileNameSmall = ['COMSOL_' axName '_' OM '_' rot suff]; % used for both loading and saving
resultsFileName = [resultsDir '/' resultsFileNameSmall '.mat'];

if loadResults == 1
    % check if file exists, load file
    if isfile(resultsFileName)
        disp(['loading ', resultsFileName])
        load(resultsFileName)
    else
        disp(['existing file named ' , resultsFileNameSmall,'.mat not found, creating new results table'])
        sparseSensorResults = table();
    end
else
    sparseSensorResults = table();
    if saveResults == 1
        disp(['Warning: This will overwrite any existing results file named ' resultsFileNameSmall '.mat'])
    end
end


%% Run simulation
cvx_setup;  % If this line does not run and you have already downloaded cvx, make sure it is in your MATLAB path

nSensors = Pars.chordElements*Pars.spanElements;

for iPar = 1:length(ParCombo)
    warning('off', 'all')

    
    % 1) Load simulated strain data generated in COMSOL
    StrainSet = process_comsol_strain(ParCombo,iPar,homedir,OM,axName,rot,tStep);
    
    
    % 2) transform simulated data with neural sensors 
    
    % find probability of firing (neural filtering and nonlinearity)
    [X_ne, G_ne, ~] = neuralTransformationOfData(StrainSet,ParCombo(iPar));
    
    for iter = 1:simIter
    
        % convert probability of firing to spiking data
        timePtsPerSpikeRep = size(X_ne,2)/(ParCombo(iPar).sampFreq / ParCombo(iPar).flapFrequency)-length(unique(G_ne)); % first spike of each condition is removed
        X = zeros(size(X_ne,1),timePtsPerSpikeRep*ParCombo(iPar).spikeReps);
        G = zeros(1,timePtsPerSpikeRep*ParCombo(iPar).spikeReps);
        
        for spRep = 1:ParCombo(iPar).spikeReps
            thisRepIdx = (spRep-1)*timePtsPerSpikeRep+1:spRep*timePtsPerSpikeRep;
            [X(:,thisRepIdx), G(thisRepIdx), ~] = convertProbFiringToSpikes(X_ne,G_ne,ParCombo(iPar));
        end
        
        % split data into train and test sets
        [X, G, XTest, GTest] = trainTestSplit(X,G,ParCombo(iPar));

        
        % 3) Dimensionality reduction on sensor data (PCA)
        
        %%% standardize training data
        Xmean = mean(X,2);
        Xstd = std(X,[],2);
        Xstd(Xstd<1e-14) = 1; % avoid divide by zero error
        XNorm = (X-Xmean)./repmat(Xstd,1,size(X,2));
        
        [w_t, Psi] = dimReductionForSspoc(XNorm, G, ParCombo(iPar));
        
        
        % 4) Optimization to find sensor locations
        
        [sensors, cutoffLim, s] = sspocOptim(w_t, Psi, ParCombo(iPar), length(unique(G)));
        [~, I_top] = sort( sum(abs(s),2),'descend');
        
        % calculate classification accuracy
        sensorsSort = I_top(1:Pars.rmodes);
        sensors10 = sensorsSort(1:10);
        
        acc10 = classAccuracyLin(X, G, XTest, GTest, sensors10);
        accAll = classAccuracyLin(X, G, XTest, GTest, 1:size(X,1));
        acc = classAccuracyLin(X, G, XTest, GTest, sensors);
        
        
        % format results for saving
        results = formatResults(sensors, acc, sensors10, acc10, accAll, ParCombo(iPar));
        sparseSensorResults = [sparseSensorResults; results]; % add new results to existing table
        
        % plot results
        if plotFlag == 1
            plotSensorLocation(sensors, acc10, ParCombo(iPar))
            figure; plotPSTHsPopulationOverlay(X,G,sensors10(1:5)',Pars.sampFreq,Pars.flapFrequency,0.1,[],1,gca);
        end
        
        disp(['parameter set ' num2str(iPar) '/' num2str(length(ParCombo)) ', simIter ' num2str(iter) '/' num2str(simIter) ' done'])
 
    end % end iteration sweep loop
    disp(['parameter set ' num2str(iPar) '/' num2str(length(ParCombo)) ' done'])

end % end parameter sweep loop


%% save results to file

if saveResults == 1
    save(resultsFileName, 'sparseSensorResults')
    disp(['saved updated results in ', resultsFileName])
end