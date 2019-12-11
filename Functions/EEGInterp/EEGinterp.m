function [Approx]=EEGinterp(RBF,c,daten,defekte_elektroden)
%Interpolates EEG data using radial basis functions 
%
%The distribution of the electrodes needs to be specified in polar coordinates
%in this example the distribution is included in the files Polarkoordinaten(2).dat
%
%The electrodes are numbered according to their appearence in the vector Polarkoordinaten.dat
%
%INPUT:
% RBF := Which radial basis function to be used, possible choices:
%        'TPS':= Thin-plate-splines / 'r^3' := Konst^3 
%        'MQ':= Multiquadric/ 'IMQ':= Inverse MQ / 'Gau':= Gaussian
% c := smoothing parameter for the radial basis functions (only needed for TPS MQ IMQ Gau)
% daten := Matrix of measurements( format:=(time X electrode position))
% defekte_elektroden := vector of electrodes to be reconstructed 
%
%OUTPUT:
% Approx := Full set of measurements with the reconstruction in places of defekte_elektroden
%
%
% Copyright (C) 2016 Janin Jäger
% 
% Justus-Liebig-Universität Gießen, Germany
% Author: janin.jaeger@math.uni-giessen.de
% Created: 4.11.2016

% Adapted to Biosemi EEG device / EEGLAB
% by Corentin Wicht / corentin.wicht@unifr.ch (06.05.2019)

% Detecting the number of channels
NumEl = daten.nbchan;

% Catch matrix size
MatSize = size(daten.data);

% May be useless?
% if length(MatSize) < 3
%    error('Insufficient Number of trials !') 
% end

% Reshaping EEGLAB data to match current format
daten.data = permute(reshape(daten.data, daten.nbchan, daten.pnts*daten.trials),[2 1]);

% Compute Distancematrix Abs
if ~isempty(daten.chaninfo.nodatchans)
    [Abs,Coord]=COS(NumEl,'refchanstruct',daten.chaninfo.nodatchans);  
else
    [Abs,Coord]=COS(NumEl);  
end

%Compute the interpolation matrix IM from Abs and the RBF
[IM,l]=RBFS(Abs,RBF,c);

%Compute if necessary polynomial part of the interpolation matrix
PolM=[];
for L=0:l-1
    for M=-L:L
PolM=[PolM,SpherHam(L,M,Coord)];
    end
end

%Add polynomial part to interpolation matrix
IMP=[IM;PolM'];

% Compute evaluation matrix EM as copy of the interoplation matrix from which
% the columns of the defekte_elektroden are erased
EM=[IM,PolM];

%Add part for the additional interpolation conditions to the interpolation matrix
O=zeros(size(PolM,2),size(PolM,2));
Q=[PolM;O];
IM=[IMP,Q];               


            
%Remove the columns of the electrodes to be reconstructed from IM, EM and daten

defekte_elektroden=sort(defekte_elektroden);                   
                    h=0;
for defekte_elektrode=defekte_elektroden

            IM(:,defekte_elektrode-h)=[];
            EM(:,defekte_elektrode-h)=[];
            IM(defekte_elektrode-h,:)=[];                    
            daten.data(:,defekte_elektrode-h)=[];
            h=h+1;
end

                  
%Produce a unit matrix with added zero rows for the cardinal solution of the problem

EP= eye(size(IM,1)-size(PolM,2),size(IM,1)-size(PolM,2));
EP=[EP;zeros(size(PolM,2),size(EP,2))];


%Solve the cardinal interpolation problem

Interpolationsmatrix=IM\EP;
                  

%Produce the approximation
        
Approx=(EM*(Interpolationsmatrix*daten.data'));

% Reshape data to match EEGLAB format
if length(MatSize) > 2
    Approx = reshape(Approx, MatSize(1), MatSize(2), MatSize(3));
end

end
