% 12/03/2018 Makoto. patch() behavior has changed after Matlab 2014. Now double(figureHandle) to obtain the figure handle.

function plot_head_region(headGrid, sumDipoleDensity, varargin)
% plot_head_region(headGrid, membershipCube, regionColor, regionOptions, showProjectedOnMrs, projectionAxis)

inputOptions = finputcheck(varargin, ...
   {'regionColor'         {'real' 'string'} [] [0.3 1 0.3]; ...
    'regionOptions'        'cell'    {} {};...
    'showProjectedOnMrs'   'boolean' [] true;
    'projectionAxis'       'real'    []  [1 2 3];...
    'projectionAlpha'      'real'    []  0.30;...
    'edgeAlpha'            'real'    []  0.60;...
    'edgeColor'            'real'    []  [0.5 0.5 0.5];...
    'edgeLinePresent'      'boolean' []  1;...
    });

if nargin < 5
    showProjectedOnMrs = true;
end;

% Obtain membershipCube. 12/04/2018 Makoto.
membershipCube = logical(sumDipoleDensity);

% Normalize dipole density so that max value == 1.
sumDipoleDensity = sumDipoleDensity/max(sumDipoleDensity(:));

position = headGrid.getPosition(membershipCube);

serializedPatchValues = [];

for p = 1:size(position,1)
    cubeSize = headGrid.spacing;
    x = position(p,1) - cubeSize * 0.5 +  [0 1 1 0 0 0;1 1 0 0 1 1;1 1 0 0 1 1;0 1 1 0 0 0]*cubeSize;
    y = position(p,2) - cubeSize * 0.5 + [0 0 1 1 0 0;0 1 1 0 0 0;0 1 1 0 1 1;0 0 1 1 1 1]*cubeSize;
    z = position(p,3) - cubeSize * 0.5 + [0 0 0 0 0 1;0 0 0 0 0 1;1 1 1 1 0 1;1 1 1 1 0 1]*cubeSize;
    for i=1:6
        
        % prevent face that are already plotted to be painted agains.
        if ~isempty(serializedPatchValues)
            potentialSerializedPatchValues = [mean(x(:,i));mean(y(:,i));mean(z(:,i))];
            difference  = serializedPatchValues(1:end-1,:) - repmat(potentialSerializedPatchValues, 1, size(serializedPatchValues,2));
            sumAbsoluteDifference = max(abs(difference));
            sameFaceId  = find(sumAbsoluteDifference < eps);
            
            % delete faces inside (between) other faces (cubes).
            if ~isempty(sameFaceId) && ishandle(serializedPatchValues(end, sameFaceId))
                delete(serializedPatchValues(end, sameFaceId));
            end;
        end;
        
        if isempty(serializedPatchValues) || ~any(sumAbsoluteDifference < eps)
            patchHandle = patch(x(:,i),y(:,i),z(:,i),'r');
                        
            %serializedPatchValues = cat(2, serializedPatchValues, [mean(x(:,i));mean(y(:,i));mean(z(:,i)); h]); % 11/30/2018 Makoto.
            serializedPatchValues = cat(2, serializedPatchValues, [mean(x(:,i));mean(y(:,i));mean(z(:,i)); double(patchHandle)]); % https://www.mathworks.com/matlabcentral/answers/356774-getting-a-numeric-graphics-handle-for-a-patch
    
            % Set visualization parameters. 12/03/2018 Makoto. Modified.
            if inputOptions.edgeLinePresent == 1
                currentLineStyle = '-';
            else
                currentLineStyle = 'none';
            end
            set(patchHandle, 'facecolor', inputOptions.regionColor,...
                             'FaceAlpha', inputOptions.projectionAlpha,...
                             'EdgeAlpha', inputOptions.edgeAlpha,...
                             'EdgeColor', inputOptions.edgeColor,...
                             'LineStyle', currentLineStyle);
            
                % This part is an experimental section to show gradation to
                % the cubes according to dipole density. It was not very
                % clean, so not adapted. 12/04/2018 Makoto.
                %
                %             % Obtain current x, y, z coordinates.
                %             xPos = mean(x(:,i));
                %             yPos = mean(y(:,i));
                %             zPos = mean(z(:,i));
                %             
                %             % Convert coordinates into indices.
                %             [xMin,xIdx] = min(abs(headGrid.xCube(1,:,1)-xPos));
                %             [yMin,yIdx] = min(abs(headGrid.yCube(:,1,1)-yPos));
                %             [zMin,zIdx] = min(abs(headGrid.zCube(1,1,:)-zPos));
                %             
                %             % Obtain the dipole density of the voxel.
                %             dataDrivenAlpha = sumDipoleDensity(yIdx, xIdx, zIdx);
                %             
                %             % This process assumes there be only one axis off the voxel at a time.
                %             if dataDrivenAlpha == 0
                %                 
                %                 if xMin ~= 0
                %                     xIdx = xIdx+1;
                %                 end
                %                 if yMin ~= 0
                %                     yIdx = yIdx+1;
                %                 end
                %                 if zMin ~= 0
                %                     zIdx = zIdx+1;
                %                 end
                %                 
                %                 % Obtain the dipole density of the voxel.
                %                 dataDrivenAlpha = sumDipoleDensity(yIdx, xIdx, zIdx);
                %                 
                %                 if dataDrivenAlpha == 0
                %                     error('Did not work')
                %                 end
                %             end
                %             
                %             set(patchHandle, 'facecolor', inputOptions.regionColor, 'FaceAlpha', dataDrivenAlpha, 'EdgeAlpha', inputOptions.edgeAlpha, 'EdgeColor', inputOptions.edgeColor, 'LineStyle', 'none');

            if nargin > 3 && ~isempty(inputOptions.regionOptions)
                set(patchHandle, inputOptions.regionOptions{:});
            end;
        end;
    end
end;

if inputOptions.showProjectedOnMrs
    plot_head_region_projected_on_mrs(headGrid, membershipCube, inputOptions.regionColor, inputOptions.regionOptions, inputOptions.projectionAxis, inputOptions.projectionAlpha);
end;