function [connec,coord,bface] = funcReadInputMesh(filepath)

fid = fopen(filepath);

if fid == -1
    error('ERROR(funcReadInputMesh): Baseline mesh .dat file not found.');
end

x   = textscan(fid,'%d',1); % Header number.
x   = textscan(fid,'%s',4); % Header strings.
ne  = textscan(fid,'%d',1); % Num elements.
np  = textscan(fid,'%d',1); % Num points.
nbf = textscan(fid,'%d',1); % Num boudary faces.

% Get connectivities.
x   = textscan(fid,'%s',1); % "Connectivities".
input = textscan(fid,'%f',(5*ne{1}));
inp = input{1};
connec = zeros(length(inp)/5,3);
for i = 1:(length(inp)/5)
    connec(i,1) = inp(5*i-3);
    connec(i,2) = inp(5*i-2);
    connec(i,3) = inp(5*i-1);
end

% Get coordinates.
x = textscan(fid,'%s',1); % "Coordinates".
input = textscan(fid,'%f',(5*np{1}));
inp = input{1};
coord = zeros(length(inp)/5,2);
for i = 1:(length(inp)/5)
    coord(i,1) = inp(5*i-3);
    coord(i,2) = inp(5*i-2);
end

% Get bounds.
x = textscan(fid,'%s%s',1); % "Boundary Faces".
input = textscan(fid,'%f',(3*nbf{1}));
inp = input{1};
bface = zeros(length(inp)/3,3);
for i = 1:(length(inp)/3)
    bface(i,1) = inp(3*i-2);
    bface(i,2) = inp(3*i-1);
    bface(i,3) = inp(3*i);
end

% % Print mesh details.
% fprintf('MESH DETAILS: \n');
% fprintf('         Number of ELEMENTS: %d\n',  ne{1});
% fprintf('           Number of POINTS: %d\n',  np{1});
% fprintf('   Number of BOUNDARY FACES: %d\n\n',nbf{1});

fclose(fid);

end