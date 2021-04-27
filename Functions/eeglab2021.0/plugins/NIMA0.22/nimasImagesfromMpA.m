% nimasImagesfromMpA() - Generates blobs or voxels inside the brain space
%                        using images from Nima's Measure Projection. Input
%                        must be single or multiple cells containing 1x3
%                        xyz coordinates of dipole locations. It applies
%                        3-D Gaussian blurr kernel to dipole locations to
%                        make them into probabilistic dipole density, and
%                        visualize it with either blobs or voxels. The
%                        preset voxels size is 8x8x8 mm (determined by
%                        headGrid.mat). This function requires several
%                        files that were stripped fromo Measure Projection
%                        Toolbox by Makoto to make it available as a
%                        stand-alone function.
%
% Usage:
%
%       nimasImagesfromMpA(dipoleLocations, varargin)
%
% Input:
%
%       dipoleLocations - {1 x n} cells, each containing [m x 3] which is
%                         m sets of EEG.dipfit.model.posxzy.
%
% Optional Input as varargin:
%
%     'voxelSize' - spatial resolution in mm. [default: 4]
%     'FWHM' - full-width half-maximum of 3-D Gaussian kernel in mm. [default: 8]
%     'numGaussTrunc' - number of sigma to truncate the Gaussian. [default: 2]
%     'blobOrVoxel' - selection of the type of visualization. [default: 'blob']
%     'uniformAlpha' - transparency of the blob/voxel. [0-1 default: 0.25]
%     'edgeAlpha' - applies to voxel plot only. [default: 0.6]
%     'edgeColor' - applies to voxel plot only. [default: 0.5 0.5 0.5]
%     'edgeLineStyle' - applies to voxel plot only. Can be 'none' to delete the voxel lines. [default: '-']
%     'mainLightIntensity' - brightness of the main light. [default: 1]     
%     'secondaryLightIntensity' - same as above. [default: 1]
%     'projectionPlanes' - 0-no MRI image or projection. 1-Sagittal, 2-Coronal, 3-Axial. [default: 1 2 3] 
%     'separateAlphaList' - overwrites 'uniformAlpha' above. This needs to 
%                           be [n x 1] where n is the number of cells in the input
%                           'dipoleLocations'.
%     'separateColorList' - this needs to be [n x 3] where n is the number of cells 
%                           in the input 'dipoleLocations' and 3 is RGB values.
%                           If not specified, 13 W3C colors will be used by default.
%                           See below for the list of the 13 colors.
%
%                   colors{1,1}  = [1 1 1];            % White (excluded)
%                   colors{2,1}  = [1 1 0];            % Yellow
%                   colors{3,1}  = [1 0 1];            % Fuchsia
%                   colors{4,1}  = [1 0 0];            % Red
%                   colors{5,1}  = [0.75  0.75  0.75]; % Silver
%                   colors{6,1}  = [0.5 0.5 0.5];      % Gray (excluded)
%                   colors{7,1}  = [0.5 0.5 0];        % Olive
%                   colors{8,1}  = [0.5 0 0.5];        % Purple
%                   colors{9,1}  = [0.5 0 0];          % Maroon
%                   colors{10,1} = [0 1 1];            % Aqua
%                   colors{11,1} = [0 1 0];            % Lime
%                   colors{12,1} = [0 0.5 0.5];        % Teal
%                   colors{13,1} = [0 0.5 0];          % Green
%                   colors{14,1} = [0 0 1];            % Blue
%                   colors{15,1} = [0 0 0.5];          % Navy
%                   colors{16,1} = [0 0 0];            % Black (excluded)

% Copyright (C) 2018, Makoto Miyakoshi and Nima Bigdely-Shamlo, SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

% History
% 12/06/2018 Makoto. Created for Christiane Glatz.

function nimasImagesfromMpA(dipoleLocations, varargin)

% Remove or add optionalInputs.
if isempty(varargin{end})
    varargin = varargin(1:end-1);
else
    vararginContents = varargin{end};
    varargin = varargin(1:end-1);
    varargin = [varargin vararginContents];
end

% Prepare inputs.
inputOptions = finputcheck(varargin, ...
    {'voxelSize'               'real'    [2 8] 4;...
     'FWHM'                    'real'    []    8;...
     'numGaussTrunc'           'real'    []    2;...
     'blobOrVoxel'             'real'    [1 2] 1;...
     'uniformAlpha'            'real'    []    0.25;...
     'separateAlphaList'       'real'    []    [];...
     'edgeAlpha'               'real'    []    0.60;...
     'edgeColor'               'real'    []    [0.5 0.5 0.5];...
     'edgeLinePresent'         'boolean' []    1;...
     'separateColorList'       'real'    []    [];...
     'mainLightIntensity'      'real'    [0 1] 1;...
     'secondaryLightIntensity' 'real'    [0 1] 1;...
     'projectionPlanes'        'integer' [0 3] [1 2 3];...
    });



%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare colors. %%%
%%%%%%%%%%%%%%%%%%%%%%%
if isempty(inputOptions.separateColorList)
    
    separateColorListFlag = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Set 12 blob colors.%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    % 16 color names officially supported by W3C specification for HTML
    colors{1,1}  = [1 1 1];            % White (excluded)
    colors{2,1}  = [1 1 0];            % Yellow
    colors{3,1}  = [1 0 1];            % Fuchsia
    colors{4,1}  = [1 0 0];            % Red
    colors{5,1}  = [0.75  0.75  0.75]; % Silver
    colors{6,1}  = [0.5 0.5 0.5];      % Gray (excluded)
    colors{7,1}  = [0.5 0.5 0];        % Olive
    colors{8,1}  = [0.5 0 0.5];        % Purple
    colors{9,1}  = [0.5 0 0];          % Maroon
    colors{10,1} = [0 1 1];            % Aqua
    colors{11,1} = [0 1 0];            % Lime
    colors{12,1} = [0 0.5 0.5];        % Teal
    colors{13,1} = [0 0.5 0];          % Green
    colors{14,1} = [0 0 1];            % Blue
    colors{15,1} = [0 0 0.5];          % Navy
    colors{16,1} = [0 0 0];            % Black (excluded)
    
    % Silver is twice brighter because the color may be used for background
    colors{5,1} = [0.875 0.875 0.875];
    
    % Exclusing:White, Black, Gray. Choosing and sorting the remaining 13 colors: Red, Blue, Green, Yellow, Fuchsia, Lime, Aqua, Maroon, Olive, Purple, Teal, Navy, and Silver.
    inputOptions.separateColorList = colors([4 13 14 2 3 11 10 9 7 8 12 15 5]);
else
    
    separateColorListFlag = 1;
    
    % Convert nx3 matrix to cell array.
    colors = cell(size(inputOptions.separateColorList, 1),1);
    for colorIdx = 1:size(inputOptions.separateColorList, 1) 
        colors{colorIdx,1} = inputOptions.separateColorList(colorIdx,:);
    end
    inputOptions.separateColorList = colors;
    
    if length(dipoleLocations) == size(inputOptions.separateColorList, 1)
        disp('User defined color will be used.')
    else
        error('The number of colors should be the same as the number of blobs to plot')
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Prepare 3-D Gaussian kernel. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct headGrid object.
headGridObj = headGrid(inputOptions.voxelSize);

% Construct parametersForProjection object for the 3-D Gaussian blurr kernel.
sigmaInGaussian = inputOptions.FWHM/2.355; % This calculates sigma in Gaussian equation.
parametersForProjectionObj = projectionParameter(sigmaInGaussian);
parametersForProjectionObj.numberOfStandardDeviationsToTruncatedGaussaian = inputOptions.numGaussTrunc;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Open the main figure. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure
plot_dipplot_with_cortex;
plotOptions   = {};
regionOptions = {};



%%%%%%%%%%%%%%%%%%%%%%%
%%% Start the loop. %%%
%%%%%%%%%%%%%%%%%%%%%%%
for numBlobIdx = 1:length(dipoleLocations)
    
    % Obtain the current dipole cluster.
    dipoleInfo.location = dipoleLocations{numBlobIdx};
    
    % Transform dipole locations into dipole density.
    [projectionMatrix, totalDipoleDenisty, dipoleDensity] = getProjectionMatrix(dipoleInfo, headGridObj, parametersForProjectionObj);
    
    % Obtain inside-brain inclusive mask.
    insideBrainInclusiveMaskIdx = find(vec(headGridObj.insideBrainCube));
    
    % Create a mask for inside-brain inclusive voxel mask on 27x23x23 tensor.
    emptyBrainToBePopulated1D = vec(zeros(headGridObj.cubeSize));
    emptyBrainToBePopulated1D(insideBrainInclusiveMaskIdx) = sum(dipoleDensity,1);
    sumDipoleDensity = reshape(emptyBrainToBePopulated1D, headGridObj.cubeSize);

    % Determine the current color index.
    if separateColorListFlag == 0;
        colorIdx = mod(numBlobIdx, length(inputOptions.separateColorList));
        if colorIdx == 0
            colorIdx = length(inputOptions.separateColorList);
        end
    else
        colorIdx = numBlobIdx;
    end
    
    % Show blob or voxel.
    switch inputOptions.blobOrVoxel
        case 1
            if isempty(inputOptions.separateAlphaList)
                surfaceOptions = {'FaceAlpha' inputOptions.uniformAlpha};
            else
                surfaceOptions = {'FaceAlpha' inputOptions.separateAlphaList(numBlobIdx)};
            end
            
            plot_head_surface(headGridObj, logical(sumDipoleDensity),...
                'surfaceColor', inputOptions.separateColorList{colorIdx},...
                'surfaceOptions', surfaceOptions,...
                'mainLightIntensity', inputOptions.mainLightIntensity,...
                'secondaryLightIntensity', inputOptions.secondaryLightIntensity,...
                'projectionAxis', inputOptions.projectionPlanes,...
                plotOptions{:});
            
        case 2
            if isempty(inputOptions.separateAlphaList)
                projectionAlpha = inputOptions.uniformAlpha;
            else
                projectionAlpha = inputOptions.separateAlphaList(numBlobIdx);
            end
            plot_head_region(headGridObj, sumDipoleDensity,...
                'regionColor',  inputOptions.separateColorList{colorIdx},...
                'projectionAlpha', projectionAlpha,...
                'edgeAlpha', inputOptions.edgeAlpha,...
                'edgeColor', inputOptions.edgeColor,...
                'edgeLinePresent', inputOptions.edgeLinePresent,...
                'regionOptions',  regionOptions,...
                'projectionAxis', inputOptions.projectionPlanes);
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Delete MRI images that are not selected. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the default order: Axial(109x91x3), Sagittal(91x109x3), Coronal(91x91x3).
defaultOrder = [109 91 3; 91 109 3; 91 91 3];

% Calculate the CDAta order.
mrImageHandles = findobj(gca, 'Tag', 'img');
CDataSize1 = size(mrImageHandles(1).CData);
CDataSize2 = size(mrImageHandles(2).CData);
CDataSize3 = size(mrImageHandles(3).CData);
CDataOrder = [find(sum(abs(bsxfun(@minus, defaultOrder, CDataSize1)),2)==0) ...
    find(sum(abs(bsxfun(@minus, defaultOrder, CDataSize2)),2)==0) ...
    find(sum(abs(bsxfun(@minus, defaultOrder, CDataSize3)),2)==0)];

% Define user input order: Sagittal, Coronal, Axial.
convertUserInputToDefaultOrder = [2 3 1];

% Calculate the image index to plot.
if inputOptions.projectionPlanes == 0
    delete(mrImageHandles([1 2 3]));
    
else
    inputOptions.projectionPlanes = nonzeros(inputOptions.projectionPlanes); % Foolproof by removing zero.
    
    imageIdxToPlot = convertUserInputToDefaultOrder(inputOptions.projectionPlanes);
    
    % Calculate the image index to delete.
    imageUridxToDelete = setdiff([1 2 3], imageIdxToPlot);
    
    % Determine the image indices to delete.
    [~, imageIdxToDelete] = intersect(CDataOrder, imageUridxToDelete);
    
    % Delete the MR images as necessary.
    delete(mrImageHandles(imageIdxToDelete));
end