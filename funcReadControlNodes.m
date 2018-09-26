function [CNs] = funcReadControlNodes(filepath)

fid = fopen(filepath);

if fid == -1
    error('ERROR(funcReadControlNodes): Control_Nodes.txt file not found.');
end

% Pull out CNs.
input = textscan(fid,'%f');
input = input{1};
CNs = zeros(length(input)/2,2);
for i = 1:(length(input)/2)
    CNs(i,1) = input(2*i-1);
    CNs(i,2) = input(2*i);
end

fclose(fid);

end