function [X, G, filtStrain] = neuralTransformationOfData(StrainSet, Pars)
% [X, G, filtStrain] = neuralTransformationOfData(StrainSet,Pars)
%
% Takes in strain data and converts it to neurally encoded strain via a
% linear-nonlinear model.
%
% Inputs:
%   StrainSet: structure with strain data
%     Fields:
%      strain_0
%      strain_1
%   Pars: structure with parameters for simulation.
%     Fields:
%      staWidth: width of STA (ms)
%      staFreq: frequency of STA (Hz)
%      staDelay: time delay before oscillations of STA begin (ms)
%      nldShift: horizontal offset (threshold shift) of NLD
%      nldGrad: slope parameter of NLD
%      normalizeVal
%      subSamp
%
% Outputs:
%   X: data matrix; nSensorLocs x nDataPts
%   G: vector indicating classes of columns in X; 1 x nDataPts
%   filtStrain: matrix of filtered strain; nSensorLocs x nDataPts

fieldsStrainSet = fields(StrainSet);
nSensorLocs = size(StrainSet.(fieldsStrainSet{1}),1);

staFunc = @(t) cos(  Pars.staFreq*(t+Pars.staDelay)) .* exp(-(t+Pars.staDelay).^2 / Pars.staWidth.^2  );
nldFunc = @(s) (  1 ./ ( 1+exp(-Pars.nldGrad.*(s-Pars.nldShift)) ) - 0.5  ) + 0.5;

staT = -19:1000/Pars.sampFreq:0;

f = staFunc(staT);
f = f-mean(f);
if Pars.staFreq < .1
    f = ones(size(f));
end

k = sqrt(1/sum(f.^2)); % k = 0.2003 for original parameters
staFilt = fliplr(k*f /0.2003 *1000/Pars.sampFreq);  

nConditions = 2;
n_conv = ( Pars.simStartup  *Pars.sampFreq*Pars.subSamp +2 -length(staT) )...  % set time points of data that will be used for convolution
    : Pars.simEnd*Pars.sampFreq*Pars.subSamp;
n_out = round((Pars.simEnd-Pars.simStartup) * Pars.sampFreq*Pars.subSamp);

convMat = zeros(nSensorLocs,n_out*nConditions);
G = zeros(1,n_out*nConditions);

for iRot = 1:nConditions     % convolve data with neural filter for each rotation rate: 0 and 1 rad/s
    sSet = StrainSet.(['strain_' num2str(iRot-1)]);
    strainConv = zeros(nSensorLocs,n_out);
    
    for locNum = 1:nSensorLocs
        strainConv(locNum,:) = conv(sSet(locNum,n_conv),staFilt,'valid');
    end
    convMat(:,(iRot-1)*n_out+1:(n_out*iRot)) = strainConv;
    G(:,(iRot-1)*n_out+1:(n_out*iRot)) = iRot;
end % end rotation loop

filtStrain = convMat;
X = nldFunc( filtStrain/Pars.normalizeVal/Pars.subSamp );  % run filtered data through nonlinear decision function
