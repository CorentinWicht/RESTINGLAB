% EEG = fitTwoDipoles(EEG, symmetryRegion, threshold)
%
% Input:
%      EEG--EEGLAB EEG structure
%      threshold--threshold value for "true" peak selection (default = 35)
%      sym_reg--symmetry region (default LRR = Large rectangular region)
%                 other admitted values:
%                 LSR = Large Square Region
%                 SRR = Small Rectangular Region
%                 SSR = Small Square Region
% Output:
%      EEG--EEGLAB EEG structure

% Copyright (C) Makoto Miyakoshi, SCCN,INC,UCSD
%               Catarina Piazza, Politecnico di Milano/Scientific Institute IRCCS Eugenio Medea
%
% Reference: 
% Piazza, Miyakoshi, Akalin-Acar, Cantiani, Reni, Bianchi, Makeig. (2016). An automated
% function for identifying EEG independent components representing bilateral source
% activity. XIV Mediterranean Conference on Medical and Biological Engineering and
% Computing 2016 IFMBE Proceedings 57:105-109.
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
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function EEG = fitTwoDipoles(EEG, symmetryRegion, threshold)

% Run Caterina's bilateral dipole finder.
sym_class = IC_sym_class(EEG,  'sym_reg', symmetryRegion, 'threshold', threshold, 'plot', 'no');
bilateralDipoleIdx = find(sym_class(2,:));

% Fit bilateral dipoles
for m = 1:length(bilateralDipoleIdx)
    firstDipolePosxyz  = EEG.dipfit.model(bilateralDipoleIdx(m)).posxyz(1,:);
    firstDipoleMomxyz  = EEG.dipfit.model(bilateralDipoleIdx(m)).momxyz(1,:);
    secondDipolePosxyz = firstDipolePosxyz.*[-1 1 1];
    secondDipoleMomxyz = firstDipoleMomxyz.*[-1 1 1];
    EEG.dipfit.model(bilateralDipoleIdx(m)).posxyz(2,:) = secondDipolePosxyz;
    EEG.dipfit.model(bilateralDipoleIdx(m)).momxyz(2,:) = secondDipoleMomxyz;
    EEG.dipfit.model(bilateralDipoleIdx(m)).select = [1 2];
    EEG = dipfit_nonlinear(EEG, 'symmetry', 'x', 'component', bilateralDipoleIdx(m));
end