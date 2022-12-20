function StrainSet = process_comsol_strain(ParCombo,iPar,homedir,OM,axis,rot,tStep)

nReps = 250;  % # times to repeat chosen wingbeat (include 2 or more cycles because one will be trimmed off by simStartup)
tKeep = 2001:1/(str2double(tStep)*ParCombo(iPar).sampFreq):2201;  % time window to keep
% Note that tKeep will need to change for custom simulations.  Should be
% selected to keep a single wingbeat

% parse parameter values to strings for filenames
E = ParCombo(iPar).E;
Eexp = floor(log10(E));
Eones = floor(E/10^Eexp);
alpha = ParCombo(iPar).alpha;
alphaOnes = floor(alpha/100);
alphaDec = alpha - alphaOnes*100;
if rem(alphaDec,10) == 0
    alphaDec = alphaDec/10;
end

% create filenames for no rotation and with rotation
fname1a = ['sorted_SingleWing_' OM '_E' num2str(Eones) 'E' num2str(Eexp)];
fname1b = ['sorted_SingleWing_' axis '_' OM '_E' num2str(Eones) 'E' num2str(Eexp)]; 
fname2 = ['_Z' num2str(alphaOnes) 'p' num2str(alphaDec)];
fname3 = ['_t' tStep '.csv'];
disp([ fname1b fname2 '_rot0p0' fname3]);

% process data without rotation
strain_0 = csvread([homedir '/simulated_wing_data/' fname1a fname2 '_rot0p0' fname3]);
strain_0 = interp1(-2:(size(strain_0,2)-3),strain_0',tKeep,'spline'); % up sampling freq; trim off startup time
strain_0 = strain_0(2:end,:);
strain_0 = repmat(strain_0',1,nReps);
StrainSet.strain_0 = strain_0;

% process data with rotation
strain_1 = csvread([homedir '/simulated_wing_data/' fname1b fname2 '_rot' rot fname3]);
strain_1 = interp1(-2:(size(strain_1,2)-3),strain_1',tKeep,'spline'); % up sampling freq; trim off startup time
strain_1 = strain_1(2:end,:);
strain_1 = repmat(strain_1',1,nReps);
StrainSet.strain_1 = strain_1;

