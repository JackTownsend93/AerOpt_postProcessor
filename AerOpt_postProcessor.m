%% ----------------------- AEROPT POST-PROCESSOR ----------------------- %%
clear; clc; close all;
fprintf('%%---------- AerOpt POST-PROCESSOR ----------%%\n');
fprintf('Enter user parameters and filepath to outputs\n')
fprintf('in the USER INPUTS section.\n\n');

%% USER INPUTS.
% Make sure case & run names are correct, with each string terminating in 
% a "/" character.
% E.g.: caseName = 'XC10_preOpt/' and runName = 'AerOpt2D_3.1_180820_1219/'
caseName   = 'aerOpt_xc10_fixedRearWing/';
 runName   = 'AerOpt2D_3.1_180830_1504/';

% Post-pro options.
nestToPlot = 1;            % Set to 1 to plot best nest, 2 for next, etc.
plotArea   = [-9 0 0 2.5]; % Zoom all plots to area of interest.
numBodies  = 5;            % May need to be bigger than expected in event 
                           % of surface jumping ( e.g.: at thin parts of 
                           % airfoils).
                           
% Boundary-plotting options (tricky).
% Area in which bodies of interest lay, make sure to exclude farfield
% boundaries.
areaOfBodies = [-9 0 0.001 2.5];
% Set start point for knn search. Try not to start this near thin sections
% where jumping is likely to occur. For XC10, place at [-8.5 1.5].
searchPoint  = [-8.5 1.5];
% Distance between knn jumps for which it should be considered to have 
% jumped to a new body.
jumpCriteria = 0.1;
                           
% Read parameters from 'AerOpt_InputParameters.txt' file.
filepath = [caseName,'Input_Data/AerOpt_InputParameters.txt'];
[Ma,Tamb,Pamb,R,gamma,Re,Low2Top,NoSnap,NoNests,NoCN,DoF,NoG,objFunc,...
 turbModel,NoSolIter,meshMove,baselineMesh] =funcReadInputParams(filepath);


%% WRITE MOST RELEVANT PARAMETERS.
fprintf('AERODYNAMIC PARAMS:\n');
fprintf('      Mach Number:    %.2f\n'   , Ma);
fprintf('     Ambient Temp:    %.0f K\n' , Tamb);
fprintf(' Ambient Pressure:    %.0f Pa\n', Pamb);
fprintf('     Reynolds Num:    %.3e\n'   , Re);
fprintf(' Turbulence Model:    %s.\n'    , turbModel);
fprintf(' Num Solver Iters:    %.0f\n'   , NoSolIter);
fprintf('\nOPTIMISATION PARAMS:\n');
fprintf('   Low:High Nests:    %.3f\n'   , Low2Top);
fprintf('    Num Snapshots:    %.0f\n'   , NoSnap);
fprintf('        Num Nests:    %.0f\n'   , NoNests);
fprintf('          Num CNs:    %.0f\n'   , NoCN);
fprintf('         Num DoFs:    %.0f\n'   , DoF);
fprintf('  Num Gens, Limit:    %.0f\n'   , NoG);
fprintf('     Obj function:    %s.\n'    , objFunc);
fprintf('    Mesh Movement:    %s.\n'    , meshMove);
fprintf('    Baseline Mesh:    %s\n\n'   , baselineMesh);

%% OFFER USER A GIF OPTION.
gifOption = input('Create GIF of mesh movement? (y/n)\n','s');
if strcmpi(gifOption, 'y')
    fprintf('\nGIF "mesh.gif" will be saved to:\n');
    fprintf('    %s%s\n',caseName,runName);
elseif strcmpi(gifOption, 'n')
    fprintf('\nNo GIF requested.\n');
else
    error('ERROR: Must enter "y" or "n".')
end

%% OFFER USER TO TRACK DEVELOPMENT OF A CONTROL NODE.
trackCNOption = input('\nEnter CN # to track evolution for (press ENTER to not track any).\n','s');
if isempty(trackCNOption) == 1
    fprintf('\nNot tracking CN evolution.\n');
elseif isnan(str2double(trackCNOption)) == 1
    error('ERROR: NaN, ensure an integer is entered.');
elseif str2double(trackCNOption) > NoCN || str2double(trackCNOption) < 0
    error('ERROR: Number lies outside of valid CN # range.');
else
    trackCNOption = str2double(trackCNOption);
    fprintf('Tracking evolution of Control Node %d',trackCNOption);
end

%% PROCESS FITNESS.
% Read from fitness file, including how many generations have actually run.
[fitness,fitnessBaseline,NoG_actual] = funcReadFitness(caseName,runName,NoNests,Ma);

% Print number of generations run so far based on fitness file.
fprintf('\nNumber of generations simulated so far:\n');
fprintf('  %d/%d\n\n',NoG_actual,NoG);

% Plot fitness.
f = 1;
figure(f); f = f+1; hold on; grid on;
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5]);
movegui(f-1,'northwest');
line([0 NoG_actual],[fitnessBaseline fitnessBaseline],'LineWidth',1,'Color','k');
xaxis = 1:NoG_actual;
for i = 1:NoNests
    plot(xaxis, fitness(:,i))
end
title('Fitness Development');
xlabel('Number of Generations'); ylabel('Fitness');
legendString = cell(NoNests+1,1);
legendString(1) = cellstr('Baseline');
for i = 1:NoNests
    legendString(i+1) = cellstr(['Nest ',num2str(i)]);
end
legend(legendString, 'Location','eastoutside');

%% READ AND PROCESS CN INFORMATION.
% Read from AerOpt_InputParameters.txt.
filepath = [caseName,'Input_Data/AerOpt_InputParameters.txt'];
[CNxrange,CNyrange] = funcReadCNs(filepath);

% Read from Control_Nodes.txt.
filepath = [caseName,'Input_Data/','Control_Nodes.txt'];
[CNs] = funcReadControlNodes(filepath);

% Find CNs that are fixed in both x and y.
j = 1; k = 1;
for i = 1:length(CNxrange)
    if CNxrange(i,1)==0 && CNxrange(i,2)==0 && CNyrange(i,1)==0 && CNyrange(i,2)==0
        % Thus fixed CNs.
        CNsFixedIdx(j) = i; % Save index.
        j = j+1;
    else
        % Thus free CNs.
        CNsFreeIdx(k) = i;  % Save index.
        k = k+1;
    end
end
% Save to their own matrices for simplicaity later.
CNsFixed = CNs(CNsFixedIdx,:);
CNsFree  = CNs(CNsFreeIdx, :);

% Generate coords for squares defining bounding box for each free CN.
CNxBoxCoords = zeros(5,length(CNsFreeIdx));
CNyBoxCoords = zeros(5,length(CNsFreeIdx));
for i = 1:length(CNsFreeIdx)
    % Creating coords for square from top left moving counter-clockwise.
    CNxBoxCoords(:,i) = [CNsFree(i,1)+CNxrange(CNsFreeIdx(i),1), CNsFree(i,1)+CNxrange(CNsFreeIdx(i),1), CNsFree(i,1)+CNxrange(CNsFreeIdx(i),2), CNsFree(i,1)+CNxrange(CNsFreeIdx(i),2), CNsFree(i,1)+CNxrange(CNsFreeIdx(i),1)];
    CNyBoxCoords(:,i) = [CNsFree(i,2)+CNyrange(CNsFreeIdx(i),1), CNsFree(i,2)+CNyrange(CNsFreeIdx(i),2), CNsFree(i,2)+CNyrange(CNsFreeIdx(i),2), CNsFree(i,2)+CNyrange(CNsFreeIdx(i),1), CNsFree(i,2)+CNyrange(CNsFreeIdx(i),1)];
end

%% SNAPSHOT ANALYSIS.

% Call function here.

% Plot Nests progression.


%% PLOT BASELINE MESH.
% Read initial baseline mesh.
filepath = [caseName,'Input_Data/',baselineMesh];
[connec,coord,bface] = funcReadInputMesh(filepath);

% Plot BASELINE MESH.
figure(f); f = f+1;
hold on; axis equal; grid on; axis(plotArea);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5]);
movegui(f-1,'northeast');
title('Baseline Mesh.'); xlabel('x-coords'); ylabel('y-coord');
triplot(connec,coord(:,1),coord(:,2),'-k','linewidth',0.01);
% Add CNs and their bounds to the plot.
scatter(CNsFixed(:,1),CNsFixed(:,2),'ko','MarkerFaceColor','g')
scatter(CNsFree(:,1),CNsFree(:,2),'rx','LineWidth',2)
for i = 1:length(CNsFreeIdx)
    plot(CNxBoxCoords(:,i),CNyBoxCoords(:,i),'r.-.')
end

% Mimic plot formatting to account for additional CN bound entry in legend.
h = zeros(4,1);
h(1) = plot(NaN,NaN,'k');
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor','g');
h(3) = plot(NaN,NaN,'rx','LineWidth',2);
h(4) = plot(NaN,NaN,'r.-.');
legend(h,'Mesh','CNs Fixed','CNs Free','CN bounds');

%% PLOT BASELINE BOUNDARY.
[bodiesForPlot] = funcOrderBoundaries(bface,coord,numBodies,areaOfBodies,searchPoint,jumpCriteria);

figure(f); f=f+1;
hold on; axis equal; grid on; axis(plotArea);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5]);
movegui(f-1,'southwest');
for i = 1:length(bodiesForPlot)
    xy = bodiesForPlot{1,i};
    plot(xy(:,1),xy(:,2),'k','LineWidth',1.5);
    xy = [];
end

%% PLOT OPTIMISED BOUNDARY.
% Read mesh of top nest at final generation (i.e.: OPTIMISED GEOMETRY).
filepath = [caseName,runName,'Mesh_Data/Geometry',num2str(nestToPlot),...
                '/Geometry',num2str(nestToPlot),'.dat'];
[connec,coord,bface] = funcReadMesh(filepath);

% Plot OPTIMISED BOUNDARY on same plot as baseline boundary.
[bodiesForPlot] = funcOrderBoundaries(bface,coord,numBodies,areaOfBodies,searchPoint,jumpCriteria);
for i = 1:length(bodiesForPlot)
    xy = bodiesForPlot{1,i};
    plot(xy(:,1),xy(:,2),'m');
    xy = [];
end

% Add CNs and their bounds to the plot.
scatter(CNsFixed(:,1),CNsFixed(:,2),'ko','MarkerFaceColor','g')
scatter(CNsFree(:,1),CNsFree(:,2),'rx','LineWidth',2)
for i = 1:length(CNsFreeIdx)
    plot(CNxBoxCoords(:,i),CNyBoxCoords(:,i),'r.-.')
end
% Mimic plot formatting to account for additional CN bound entry in legend.
h = zeros(4,1);
h(1) = plot(NaN,NaN,'k','LineWidth',1.5);
h(2) = plot(NaN,NaN,'m');
h(3) = plot(NaN,NaN,'ko','MarkerFaceColor','g');
h(4) = plot(NaN,NaN,'rx','LineWidth',2);
h(5) = plot(NaN,NaN,'r.-.');
legend(h,'Baseline Geometry','Optimised Goemetry','Fixed CNs','Free CNs','CN bounds');
title('Baseline and Optimised Geometries.');


%% PLOT OPTIMISED MESH WITH CNs.
figure(f); f=f+1;
hold on; axis equal; grid on; axis(plotArea);
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0, 0.5, 0.5]);
movegui(f-1,'southeast');
triplot(connec,coord(:,1),coord(:,2),'-k','linewidth',0.01);
title('Optimised Mesh.'); xlabel('x-coords'); ylabel('y-coord');
% Add CNs and their bounds to the plot.
scatter(CNsFixed(:,1),CNsFixed(:,2),'ko','MarkerFaceColor','g')
scatter(CNsFree(:,1),CNsFree(:,2),'rx','LineWidth',2)
for i = 1:length(CNsFreeIdx)
    plot(CNxBoxCoords(:,i),CNyBoxCoords(:,i),'r.-.')
end
% Mimic plot formatting to account for additional CN bound entry in legend.
h = zeros(4,1);
h(1) = plot(NaN,NaN,'k');
h(2) = plot(NaN,NaN,'ko','MarkerFaceColor','g');
h(3) = plot(NaN,NaN,'rx','LineWidth',2);
h(4) = plot(NaN,NaN,'r.-.');
legend(h,'Mesh','CNs Fixed','CNs Free','CN bounds');

%% WRITE AND SAVE MESH GIF (IF USER REQUESTS IT).
if strcmpi(gifOption, 'y')
    for i = 1:NoG_actual
        fileID = fopen([caseName,runName,'Output_Data/Geometry',num2str(i),'.dat']);
        x   = textscan(fileID,'%d',1);
        x   = textscan(fileID,'%s',6);
        ne  = textscan(fileID,'%d',1);
        np  = textscan(fileID,'%d',1);
        nbf = textscan(fileID,'%d',1);
        nbf = nbf{1};
        x   = textscan(fileID,'%s',1);
        input = textscan(fileID,'%f',(4*ne{1}));
        inp = input{1};
        connec = zeros(length(inp)/4,3);
        for j = 1:(length(inp)/4)
            connec(j,1) = inp(4*j-2);
            connec(j,2) = inp(4*j-1);
            connec(j,3) = inp(4*j);
        end
        x = textscan(fileID,'%s',1);
        input = textscan(fileID,'%f', (3*np{1}));
        inp = input{1};
        coord = zeros(length(inp)/3,2);
        for j = 1:(length(inp)/3)
            coord(j,1) = inp(3*j-1);
            coord(j,2) = inp(3*j);
        end
    %Plot.
    fig = figure;
    gifFileName = ([caseName,runName,'mesh.gif']);
    triplot(connec,coord(:,1),coord(:,2),'-k','linewidth',0.01); hold on;
    xlabel('x-coords'); ylabel('y-coords');
    title(['Gen#: ',num2str(i)]);
    daspect([1 1 1]);
    set(gcf,'units','points','position',[0,100,1500,500]);
    axis(plotArea)
    drawnow
    
    % Capture plot as image.
    frame = getframe(fig);
    image = frame2im(frame);
    [imageIndex,colourMap] = rgb2ind(image,256);
    
    % Write to end of gif file.
    if i == 1
        % Create on first iteration.
        imwrite(imageIndex,colourMap,gifFileName,'gif','Loopcount',inf);
    else
        % Append on subsequent iterations.
        imwrite(imageIndex,colourMap,gifFileName,'gif','WriteMode','append');
    end
    close(fig);
    fclose(fileID);
    end
end
