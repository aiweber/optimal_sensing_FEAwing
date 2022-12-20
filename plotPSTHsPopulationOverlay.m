function h = plotPSTHsPopulationOverlay(X,G,locs,sampFreq,flapFreq,timeRes,allColors,stdBar,axH)
%
% Plots PSTHs for all classes for sensor at a particular location.
%
% Inputs:
%  X: data matrix of spike times [nSensorLocs x timePoints]
%  G: vector of classes [1 x timePoints]
%  locs: sensor locations to plot PSTHs for
%  sampFreq: sampling frequency (points/sec)
%  flapFreq: flapping frequency (flaps/sec)
%  timeRes: time resolution for PSTH, ms (default: 0.5 ms)
%  allColors: matrix of color vectors (optional)
%  stdBar: flag to indicate whether standard deviation bar should be
%  plotted
%  axH: axis handle for plotting (optional)

if ~exist('timeRes','var') || isempty(timeRes)
    timeRes = 0.5;
end
if ~exist('allColors','var') || isempty(allColors)
    allColors = [89 57 105;
        54 104 153;
        69 151 164;
        167 195 169;
        232 200 21;
        229 82 114;
        205 94 53;
        133 113 81]/255;
end
if ~exist('stdBar','var')
    stdBar = 1;
end
if ~exist('axH','var') || isempty(axH)
    figure; hold on;
else
    axes(axH); hold on;
end

if isrow(locs)
    locs = locs';
end

nClasses = length(unique(G));

Xlocs = X(locs,:)/sampFreq*1000;  % convert spike times at this loc to ms
edges = min(min(Xlocs))-timeRes*3/2:timeRes:max(max(Xlocs))+timeRes*3/2;
centers = edges(1:end-1)+timeRes/2;
per = 1/flapFreq*1000;
yIdx = 0;

if exist('axH','var') && ~isempty(axH)
    axes(axH)
end
hold on;

for locIdx = 1:length(locs)
    colorIdx = rem(locIdx,size(allColors,1));
    if colorIdx == 0
        colorIdx = size(allColors,1);
    end
    for classNum = 1:nClasses
        %     firstTrialIdx = find(G==classNum,1,'first');
        %     nTrialsThisClass = sum(G==classNum);
        thisClassIdx = find(G==classNum);
        if nClasses == 1
            tintFacts = 0;
        else
            tintFacts = 0.4/(nClasses-1);
        end
   
        
        thisClassData = Xlocs(locIdx,thisClassIdx);
        n = histcounts(thisClassData,edges);
        % account for wraparound at times 0 and per
        idxZero = find(edges(1:end-1)<0 & edges(2:end)>0);
        idxMax = find(edges(1:end-1)<per & edges(2:end)>per);
        if ~isempty(idxZero) && ~isempty(idxMax)
            n(idxZero) = n(idxZero) + n(idxMax);
            n(idxMax) = n(idxZero);
        end
        if classNum == 1  % solid
        h.(['h' num2str(locIdx) num2str(classNum)]) = ...
            patch([centers centers(end:-1:1)],[n/max(n)*.8 zeros(size(n))]+yIdx,...
            allColors(colorIdx,:),...
            'edgecolor','none','facealpha',.6);
        elseif classNum == 2  % outline only
            idxPlotOutline = n~=0;
            idxPlotOutline(logical([abs(diff(n~=0)) 0])) = 1;
            idxPlotOutline(logical([0 abs(diff(n~=0))])) = 1;
            centersTemp = centers(idxPlotOutline);
         h.(['h' num2str(locIdx) num2str(classNum)]) = ...
            patch([centersTemp centersTemp(end:-1:1)],[n(idxPlotOutline)/max(n)*.8 zeros(size(n(idxPlotOutline)))]+yIdx,...
            'w',...
            'edgecolor',allColors(colorIdx,:),'facealpha',0);
        else
            disp('This function only handles 2 classes')
            return
        end
        
        if stdBar
            [m, s] = circularStats(thisClassData,per);
            plot(m,yIdx+.9,'.','markersize',15,'color',allColors(colorIdx,:));
            plot(m*[1 1] + s*[-1 1], (yIdx+.9)*[1 1],'linewidth',1.5,'color',allColors(colorIdx,:));
            plot(m*[1 1]+per + s*[-1 1], (yIdx+.9)*[1 1],'linewidth',1.5,'color',allColors(colorIdx,:));
            plot(m*[1 1]-per + s*[-1 1], (yIdx+.9)*[1 1],'linewidth',1.5,'color',allColors(colorIdx,:));
        end
    end
    yIdx = yIdx+1;
  
    if classNum<nClasses
        %         plot([0 per],[yIdx-0.05 yIdx-0.05],'k--')
    end
    
end

set(gca,'ytick',.9:yIdx,'yticklabel',locs)
box on
xlim([0 per])
title('PSTHs for all classes at chosen sensor locations')
ylabel('sensor location')
xlabel('time (ms)')



