function [CNxrange,CNyrange] = funcReadCNs(filepath)
% Note: relies on a consistent format for the input file, with data always
% on the same line.

% "filepath" directs to "AerOpt_InputParamters.txt" file.
fid = fopen(filepath,'r');

if fid == -1
    error('ERROR(funcReadCNs): AerOpt_InputParameters.txt file not found.');
end

%% CN range values.
% Skip first 10 lines.
for i = 1:10
    line = fgetl(fid);
end

% Pick out line 11 and 12.
xrangeLine = fgetl(fid);
yrangeLine = fgetl(fid);

% Pick out number values from string.
xrangeSplit = strsplit(xrangeLine);
yrangeSplit = strsplit(yrangeLine);
xrange = zeros(1,length(xrangeSplit)-2);
yrange = zeros(1,length(yrangeSplit)-2);
for i = 1:length(xrangeSplit)-2
    xrange(i) = str2double(xrangeSplit{i+2});
    yrange(i) = str2double(yrangeSplit{i+2});
end

% Organise x and y ranges into matrix of length = num of CNs.
CNxrange = zeros(length(xrange)/2,2);
CNxrange(:,1) = xrange(1:length(xrange)/2);
CNxrange(:,2) = xrange((length(xrange)/2+1):length(xrange));
CNyrange = zeros(length(yrange)/2,2);
CNyrange(:,1) = yrange(1:length(yrange)/2);
CNyrange(:,2) = yrange((length(yrange)/2+1):length(yrange));

fclose(fid);

end