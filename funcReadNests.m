function [topNest_xCoords,topNest_yCoords] = funcReadNests(caseName,runName,trackCN,Ma,NoCN,NoNests,NoG_actual,CNs)
%% Function to read Snapshots file.
% caseName = 'aerOpt_xc10_fixedRearWing/';
% runName  = 'AerOpt2D_3.1_180830_1504/';
% Ma       = 0.15;
% NoCN     = 20;
% NoNests  = 20;
% trackCNOption = 15;
% NoG_actual = 58;

fid = fopen([caseName,runName,'Nests',num2str(Ma,1),'.txt']);

if fid == -1
    error('ERROR(funcReadSnapshots): Nests%0.1f.txt file not found.',Ma);
end

% Initiate and pre-allocate.
count     = 0;
snapshots = zeros(NoG_actual*4*NoCN,NoNests);


for i = 1:NoG_actual
    line = fgetl(fid); % Read header line.
    for j = 1:4*NoCN
        line = fgetl(fid);                       % Get line.
        line = strsplit(line);                   % Split by spaces.
        line = line(~cellfun('isempty',line));   % Remove empty cells.
        snapshots(count+j,:) = str2double(line); % Store line in matrix.    
    end
    count = count+j;
end

% Pick out lines for user-defined CN.
% Tracking user-defined CN for the TOP NEST.
count = 0;
topNest_xCoords = zeros(NoG_actual,1);
topNest_yCoords = zeros(NoG_actual,1);

for i = 1:NoG_actual
    topNest_xCoords(i,1) = snapshots(count+trackCN     ,1);
    topNest_yCoords(i,1) = snapshots(count+trackCN+NoCN,1);
    count = count+4*NoCN;
end

% Place at CN position.
topNest_xCoords = topNest_xCoords+CNs(trackCN,1);
topNest_yCoords = topNest_yCoords+CNs(trackCN,2);

fclose(fid);

end
