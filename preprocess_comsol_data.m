%% Convert raw data files into clean data with header
%
% See https://doi.org/10.5061/dryad.fxpnvx0wq for example data files.
%
% Note that you will need to change lines 33-39 (below comment "construct
% string for raw data filename") if you create additional COMSOL data and
% use a different file naming convention.

% change the following line to point to directory with COMSOL data file
pathToComsolRawData = [pwd '/'];

% change the following line to point to directory with Github repository
pathToGithubRepo = [pwd '/'];
pathToSortedData = [pathToGithubRepo 'simulated_wing_data/'];

oms = {'0OM' '2OM'};  % order of magnitude of stiffness gradient ('0OM' or '2OM')
stiffs = { '3.0E9'};  % average stiffness, in Pa
damps = {'1.0'};      % damping ratio 
rots = {'0' '1'};     % rotation rate, in rad/s
axName = 'Yaw';       % axis of rotation
tStep = '2e-4';       % time step of COMSOL simultion, in s

for omNum = 1:length(oms)
    om = oms{omNum};
    for stiffNum = 1:length(stiffs)
        stiff = stiffs{stiffNum};
        for dampNum = 1:length(damps)
            damp = damps{dampNum};
            for rotNum = 1:length(rots)
                rot = [rots{rotNum} 'p0'];
                
                % construct string for raw data filename
                raw_file = [pathToComsolRawData ...
                    'Strain_' om ...
                    '_E' stiff ...
                    '_Zeta' damp(1) '.' damp(3) ...
                    '_' axName ...
                    '_' rot ...
                    '_t' tStep '.csv'];
                
                % construct string for sorted data filename
                if strcmp(rot,'0p0')
                    clean_file = [pathToSortedData ...
                        'sorted_SingleWing_' ...
                        om ...
                        '_E' stiff(1) stiff(4:end) ...
                        '_Z' damp(1) 'p' damp(3) ...
                        '_rot' rot ...
                        '_t2e-4.csv'];
                else
                    clean_file = [pathToSortedData ...
                        'sorted_SingleWing_' axName '_' ...
                        om ...
                        '_E' stiff(1) stiff(4:end) ...
                        '_Z' damp(1) 'p' damp(3) ...
                        '_rot' rot ...
                        '_t' tStep '.csv'];
                end
                
                
                data = csvread(raw_file,9,0);
                
                % sort data
                data(:,1:3) = round(data(:,1:3),6);
                data =sortrows(data,[1 2]);
                data = flipud(data);
                
                % save data
                dlmwrite(clean_file,data,'-append','precision',16);
            end
        end
    end
end
        
