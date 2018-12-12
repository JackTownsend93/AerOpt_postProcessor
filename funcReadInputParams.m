function [Ma,Tamb,Pamb,R,gamma,Re,Low2Top,NoSnap,NoNests,NoCN,DoF,NoG,...
    objFunc,turbModel,NoSolIter,meshMove,baselineMesh] = funcReadInputParams(filepath)

% "filepath" directs to "AerOpt_InputParamters.txt" file.
fid = fopen(filepath,'r');

if fid == -1
    error('ERROR(funcReadInputParams): AerOpt_InputParameters.txt file not found.');
end

% Search for parameter name then extract the relevant number or string
% following that regexp. Note format in InputParamaters file should be:
%
%   IV%ParamName = <number or string>   ! Description of parameter.
%
% Note the spaces between the equals are important as the third space
% delimited item is taken.

% Mach number (Ma).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'Ma');
    if lineIndex
        Ma = strsplit(line);
        Ma = str2double(Ma{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Ambient temperature (Tamb).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'Tamb');
    if lineIndex
        Tamb = strsplit(line);
        Tamb = str2double(Tamb{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Ambient pressure (Pamb).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'Pamb');
    if lineIndex
        Pamb = strsplit(line);
        Pamb = str2double(Pamb{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Specific gas constant (R).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'R');
    if lineIndex
        R = strsplit(line);
        R = str2double(R{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Specific heat ratio (gamma).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'gamma');
    if lineIndex
        gamma = strsplit(line);
        gamma = str2double(gamma{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Reynolds number (Re).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'Re');
    if lineIndex
        Re = strsplit(line);
        Re = str2double(Re{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Fraction of low to top cuckoo nests (Low2Top).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'Low2Top');
    if lineIndex
        Low2Top = strsplit(line);
        Low2Top = str2double(Low2Top{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Number of snapshots (NoSnap).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'NoSnap');
    if lineIndex
        NoSnap = strsplit(line);
        NoSnap = str2double(NoSnap{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Number of Nests (NoNests).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'NoNests');
    if lineIndex
        NoNests = strsplit(line);
        NoNests = str2double(NoNests{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Number of CNs (NoCN).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'NoCN');
    if lineIndex
        NoCN = strsplit(line);
        NoCN = str2double(NoCN{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Degrees of Freedom (DoF).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'DoF');
    if lineIndex
        DoF = strsplit(line);
        DoF = str2double(DoF{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Number of Generations (NoG).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'NoG');
    if lineIndex
        NoG = strsplit(line);
        NoG = str2double(NoG{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Objective function type (objFunc).
% Set possible options.
objFuncOptions = {'Lift/Drag','Distortion','Maximum Lift','Minimum Drag','Maximum Downforce','Minimum Lift'};
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'ObjectiveFunction');
    if lineIndex
        objFunc = strsplit(line);
        objFunc = str2double(objFunc{3});
        objFunc = objFuncOptions{objFunc};
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Turbulence model (turbModel).
% Set possible options.
turbModelOptions = {'Spalart-Allmaras','K-epsilon','SST'};
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'turbulencemodel');
    if lineIndex
        turbModel = strsplit(line);
        turbModel = str2double(turbModel{3});
        turbModel = turbModelOptions{turbModel};
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Number of solver iterations (NoSolIter).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'NoIter');
    if lineIndex
        NoSolIter = strsplit(line);
        NoSolIter = str2double(NoSolIter{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Mesh movement parameter (meshMove).
% Set possible options.
meshMoveOptions = {'Linear with Smoothing','FDGD with Smoothing','RBF','FDGD, no Smoothing'};
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'MeshMovement');
    if lineIndex
        meshMove = strsplit(line);
        meshMove = str2double(meshMove{3});
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

% Baseline mesh name (baselineMesh).
while ~feof(fid)
    line = fgetl(fid);
    lineIndex = regexp(line,'MeshFileName');
    if lineIndex
        baselineMesh = strsplit(line,'''');
        baselineMesh = [baselineMesh{2},'.dat'];
        break
    else
        % Not correct line, keep searching.
    end
end
frewind(fid);

fclose(fid);

end