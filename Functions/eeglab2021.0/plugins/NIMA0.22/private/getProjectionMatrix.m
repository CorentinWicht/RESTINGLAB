function [projectionMatrix totalDipoleDenisty gaussianWeightMatrix]= getProjectionMatrix(dipole, headGrid, projectionParameter, regionOfInterestCube)
% gaussianWeightMatrix is dipoles x (requested) grid points. It contains dipole density
% at each grid point for each dipole.

if nargin < 4 || isempty(regionOfInterestCube)
    regionOfInterestCube = headGrid.insideBrainCube;
end;

if ischar(regionOfInterestCube) && strcmpi(regionOfInterestCube, 'all')
    regionOfInterestCube = true(headGrid.cubeSize);
end;

standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo = projectionParameter.standardDeviationOfEstimatedDipoleLocation ^ 2;

% projection matrix is number of dipoles x number of grid point inside brain volume
numberOfPointsInTheResgionOfInterest = sum(regionOfInterestCube(:));

projectionMatrix = zeros(size(dipole.location , 1), numberOfPointsInTheResgionOfInterest);
totalDipoleDenisty = zeros(1, numberOfPointsInTheResgionOfInterest);
gaussianWeightMatrix = zeros(size(dipole.location , 1), numberOfPointsInTheResgionOfInterest);
distanceFromDipoleToGridLocationMatrix = zeros(size(dipole.location , 1), numberOfPointsInTheResgionOfInterest);


% a N x 3 matrix (N is the number of grid points inside brain volume
gridPosition = [headGrid.xCube(regionOfInterestCube) headGrid.yCube(regionOfInterestCube) headGrid.zCube(regionOfInterestCube);];

if projectionParameter.normalizeInBrainDipoleDenisty
    %dipoleInBrainDensityNormalizationFactor = pr.meanProjection.calculateDipoleInBrainDensityNormalizationFactor(dipole, headGrid, projectionParameter);
    dipoleInBrainDensityNormalizationFactor = meanProjection.calculateDipoleInBrainDensityNormalizationFactor(dipole, headGrid, projectionParameter);
end;

for dipoleNumber = 1:size(dipole.location , 1)
    % first place distance in the array
    distanceFromDipoleToGridLocationMatrix(dipoleNumber,:) = sum( (gridPosition - repmat(dipole.location(dipoleNumber,:), size(gridPosition,1), 1)) .^2,2 ) .^ 0.5;
    
    normalizationFactor = 1 / (projectionParameter.standardDeviationOfEstimatedDipoleLocation ^3 * sqrt(8 * (pi^3))) ;
    gaussianWeightMatrix(dipoleNumber,:) = normalizationFactor * exp(-distanceFromDipoleToGridLocationMatrix(dipoleNumber,:).^2 / (2 * standardDeviationOfEstimatedErrorInDipoleLocationPowerTwo));
    
    % truncate the dipole denisty Gaussian at ~3 standard deviation
    gaussianWeightMatrix(dipoleNumber, distanceFromDipoleToGridLocationMatrix(dipoleNumber,:) > (projectionParameter.numberOfStandardDeviationsToTruncatedGaussaian * projectionParameter.standardDeviationOfEstimatedDipoleLocation) ) = 0;
    
    % normalize the dipole in-brain denisty (make it sum up to one)
    if projectionParameter.normalizeInBrainDipoleDenisty
        gaussianWeightMatrix(dipoleNumber,:) = gaussianWeightMatrix(dipoleNumber,:) * dipoleInBrainDensityNormalizationFactor(dipoleNumber);
    end;
end;

% normalize gaussian weights to have the sum of 1 at each grid location
for gridId = 1:size(gaussianWeightMatrix, 2)
    totalDipoleDenisty(gridId) = sum(gaussianWeightMatrix(:, gridId));
    if totalDipoleDenisty(gridId) > 0
        projectionMatrix(:, gridId) = gaussianWeightMatrix(:, gridId) / totalDipoleDenisty(gridId);
    end;
end;
end