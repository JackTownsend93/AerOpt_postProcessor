function [bodiesForPlot] = funcOrderBoundaries(bface,coord,numBodies,areaOfBodies,searchPoint,jumpCriteria)
% Sorts boundary coords into order and plots distinct bodies.

%% USER INPUTS
% numBodies:    Number of distint bodies.
% areaOfBodies: User defines a box in which the bodies sit that excludes 
%               farfield points.
xmin = areaOfBodies(1);
xmax = areaOfBodies(2);
ymin = areaOfBodies(3);
ymax = areaOfBodies(4);
% plotArea:     Set axis for plot.
% searchPoint:  Starting point for knnsearch. May need to be strategically 
%               placed depending on the geometry.
% jumpCriteria: Criteria for how big a gap between points is to define a 
%               body jump.

%% SELECTING ONLY BOUNDARY FACE COORDS.
boundCoords = zeros(2*length(bface),2); %Pre-allocate.
for i = 1:length(bface)
    boundCoords(i*2-1,:) = coord(bface(i,1),1:2);
    boundCoords(i*2,:)   = coord(bface(i,2),1:2);
end
%plot(boundCoords(:,1),boundCoords(:,2))

%% ELIMINATE FARFIELD POINTS.
%Find row index of points lying within user-defined x limits.
keepIndex_X = find(boundCoords(:,1)<xmax & boundCoords(:,1)>xmin);
x_slice  = boundCoords(keepIndex_X,:); %Only keep these.
%Find row index of points lying within user-defined y limits.
keepIndex_Y = find(x_slice(:,2)<ymax & x_slice(:,2)>ymin);
xy_slice = x_slice(keepIndex_Y,:);

bodyFaces = xy_slice; %Use only this area.

% Eliminate rows of repeated coords.
[C,ia,ic] = unique(bodyFaces(:,1:2),'rows');
bodyFaces = bodyFaces(ia,:);

%% ORDER BODIES BY NEAREST NEIGHBOUR.
indexJump = zeros(numBodies,1);
numPoints = length(bodyFaces);
bodyFacesOrdered = zeros(size(bodyFaces));
for i = 1:numPoints
    [idx,D] = knnsearch(bodyFaces, searchPoint);
    bodyFacesOrdered(i,:) = bodyFaces(idx,:); % Save to ordered matrix.
    searchPoint = bodyFaces(idx,:); % Reset search point.
    bodyFaces(idx,:) = []; % Delete row to avoid recount.
end
bodyFaces = bodyFacesOrdered; clear bodyFacesOrdered;

%% SPLIT INTO DISTINCT BODIES BY IDENTIFYING LARGE JUMPS.
j = 1; idxJump = 1;
for i = 1:numPoints-1
    diffX = (bodyFaces(i,1))-bodyFaces(i+1,1);
    diffY = (bodyFaces(i,2))-bodyFaces(i+1,2);
    absDist = sqrt(diffX^2 + diffY^2);
    if absDist > jumpCriteria
        idxJump(j) = i; % Index for a body ending.
        j = j + 1;
    end
end

%% PLOT THE BOUNDARIES OF EACH BODY.
plotTracker = 1;
bodiesForPlot{1,numBodies} = [];
for i = 1:length(idxJump)
    x = bodyFaces(plotTracker:idxJump(i),1);
    y = bodyFaces(plotTracker:idxJump(i),2);
    % Add to cell.
    bodiesForPlot{1,i} = [x,y];
    plotTracker = idxJump(i)+1;
    % USE "boundary" command to seal each body rather than plot.
    % Varies too much between cases. Gives misleading plots sometimes.
    % k = boundary(x,y,0.9);
    % plot(x(k),y(k))
    % plot(x,y)
end
% Add final body to cell.
x = bodyFaces(plotTracker:numPoints,1);
y = bodyFaces(plotTracker:numPoints,2);
bodiesForPlot{1,numBodies} = [x,y];

% Plot final body.
% plot(x,y)
% axis equal; grid on;
% axis(plotArea)

end
