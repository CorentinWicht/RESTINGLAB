
function [COS,Koordinaten] = COS(a,varargin)
  
% Computes scalarproduct of the electrode positions
% INPUT (mandatory):
% a := number of electrodes here 32 or 64 is possible
%
% The distribution of the electrodes is given in the files 
%Polarkoordinaten.dat and polarkoordinaten2.dat

% INPUT (optional):
% RefChanStruct = ;


% OUTPUT: 
% COS := Matrix of scalar products of electrode positions used as distance matrix
% Koordinaten := Vector of coordinates of the electrode in Euclidean form
%
% Justus-Liebig-Universität Gießen, Germany
% Author: Janin Jäger <janin.jaeger@math.uni-giessen.de>
% Created: 2014-08-12
% Adapted to Biosemi 64-electrodes EEG device / EEGLAB
% by Corentin Wicht / corentin.wicht@unifr.ch (06.05.2019)
RefChanStruct = [];
if nargin > 1
    for i = 1:2:length(varargin)
        Param = varargin{i};
        Value = varargin{i+1};
        if ~ischar(Param)
          error('Flag arguments must be strings')
        end
        Param = lower(Param);
        switch Param
            case 'refchanstruct'
                RefChanStruct  = Value;
            otherwise
                display (['Unknown parameter setting: ' Param])
        end
    end
end

% Retrieve Channel Reference information (if defined)
%if a~=64 || a~=128
if ~isempty(RefChanStruct)
    RefElect = RefChanStruct.urchan;
end
    % Increment to match the number of channels (e.g. 63 + REF)
%     a = a+1;
% end

% Retrieve function path
Path = mfilename('fullpath');
st = dbstack;
namestr = st.name;
Path = erase(Path,namestr);
Path = [Path 'Coordinates\'];

% Retrieve all possible Polar Coordinates files
FileList = dir([Path  '**/*' '.dat']);

if ~isempty(FileList)
    CurrentFile = FileList(contains({FileList.name},num2str(a)));
    if isempty(CurrentFile)
        error('No Polar Coordinates file defined in the EEGinterp folder for %d electrodes!',a)
    end
else   
    error('No Polar Coordinates file defined in the EEGinterp folder!')
end

Polarkoordinaten=dlmread([CurrentFile.folder '\' CurrentFile.name]);

% If REF electrode, remove Polar Coordinate
if ~isempty(RefChanStruct)
    Polarkoordinaten(RefElect,:) = [];
end

%Compute Euclidean coordinates from polar coordinates

x=cos(Polarkoordinaten(:,2)*pi/180).*sin(Polarkoordinaten(:,1)*pi/180);
y=sin(Polarkoordinaten(:,2)*pi/180).*sin(Polarkoordinaten(:,1)*pi/180);
z=cos(Polarkoordinaten(:,1)*pi/180);

%Form a matrix with the Euclidean coordinates
Koordinaten=[x,y,z];

%Compute the matrix of scalar products of the electrode positions
COS=Koordinaten*Koordinaten.';
