function [Ma,Tamb,Pamb,R,gamma,Re,Low2Top,NoSnap,NoNests,NoCN,DoF,NoG,...
    objFunc,turbModel,NoSolIter,meshMove,baselineMesh] = funcReadInputParams(filepath)

% NOTE: relies on a consistent format for the input file, with data always
% on the same line.

% "filepath" directs to "AerOpt_InputParamters.txt" file.
fid = fopen(filepath,'r');

if fid == -1
    error('ERROR(funcReadInputParams): AerOpt_InputParameters.txt file not found.');
end

% Skip some lines.
for i = 1:2
    line = fgetl(fid);
end

% Ambient parameters.
Ma    = fgetl(fid);
Ma    = strsplit(Ma);
Ma    = str2double(Ma{3});

Tamb  = fgetl(fid);
Tamb  = strsplit(Tamb);
Tamb  = str2double(Tamb{3});

Pamb  = fgetl(fid);
Pamb  = strsplit(Pamb);
Pamb  = str2double(Pamb{3});

R     = fgetl(fid);
R     = strsplit(R);
R     = str2double(R{3});

gamma = fgetl(fid);
gamma = strsplit(gamma);
gamma = str2double(gamma{3});

Re    = fgetl(fid);
Re    = strsplit(Re);
Re    = str2double(Re{3});

% Skip some lines.
for i = 1:10
    line = fgetl(fid);
end

% Optimisation parameters.
Low2Top = fgetl(fid);
Low2Top = strsplit(Low2Top);
Low2Top = str2double(Low2Top{3});

NoSnap  = fgetl(fid);
NoSnap  = strsplit(NoSnap);
NoSnap  = str2double(NoSnap{3});

NoNests = fgetl(fid);
NoNests = strsplit(NoNests);
NoNests = str2double(NoNests{3});

NoCN    = fgetl(fid);
NoCN    = strsplit(NoCN);
NoCN    = str2double(NoCN{3});

line    = fgetl(fid);

DoF     = fgetl(fid);
DoF     = strsplit(DoF);
DoF     = str2double(DoF{3});

NoG     = fgetl(fid);
NoG     = strsplit(NoG);
NoG     = str2double(NoG{3});

line    = fgetl(fid);

objFunc = fgetl(fid);
objFunc = strsplit(objFunc);
objFunc = str2double(objFunc{3});
objFuncOptions = {'Lift/Drag','Distortion','Maximum Lift','Minimum Drag','Maximum Downforce','Minimum Lift'};
objFunc = objFuncOptions{objFunc};

% Skip some lines.
for i = 1:7
    line = fgetl(fid);
end

% Solver parameters.
turbModel = fgetl(fid);
turbModel = strsplit(turbModel);
turbModel = str2double(turbModel{3});
turbModelOptions = {'Spalart-Allmaras','K-epsilon','SST'};
turbModel = turbModelOptions{turbModel};

NoSolIter = fgetl(fid);
NoSolIter = strsplit(NoSolIter);
NoSolIter = str2double(NoSolIter{3});

% Skip some lines.
for i = 1:6
    line = fgetl(fid);
end

% Mesh deformation parameters.
meshMove = fgetl(fid);
meshMove = strsplit(meshMove);
meshMove = str2double(meshMove{3});
meshMoveOptions = {'Linear with Smoothing','FDGD with Smoothing','RBF','FDGD, no Smoothing'};
meshMove = meshMoveOptions{meshMove};

% Skip some lines.
for i = 1:12
    line = fgetl(fid);
end

% User input strings.
baselineMesh = fgetl(fid);
baselineMesh = strsplit(baselineMesh,'''');
baselineMesh = [baselineMesh{2},'.dat'];

fclose(fid);

end