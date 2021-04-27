% std_readdipoles() - returns cluster dipoles. I'm surprised EEGLAB did
%                     not have this simple thing. No wonder I needed to go
%                     though such a hard way to develop std_envtopo(). This
%                     function is a copy from std_dipplot() line 121-157.
% Usage:
%       clusterDipoles = std_readdipoles(STUDY, ALLEEG, [3 5 11:13])
%
% Input:
%       STUDY - EEGLAB STUDY.
%       ALLEEG - EEGLAB ALLEEG.
%       clusterIdxList - row vector of cluter indices
%
% Output:
%       clusterDipoles - cells containing the list of cluster dipoles (no
%                        centroid included). 

% Copyright (C) 2018, Makoto Miyakoshi, Hilit Serby, Arnaud Delorme, Scott Makeig. SCCN, INC, UCSD.
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
% 12/06/2018 Makoto. Disabled the part where controid is computed. Probably handing dual dipoles fails in some cases? 
% 12/04/2018 Makoto. Created for Christiane Glatz.

function clusterDipoles = std_readdipoles(STUDY, ALLEEG, clusterIdxList)

clusterDipoles = cell(1, length(clusterIdxList));

for clsIdxIdx = 1:length(clusterIdxList) % For each cluster requested
    max_r = 0;
    clear cluster_dip_models;
    len = length(STUDY.cluster(clusterIdxList(clsIdxIdx)).comps);
    ndip = 0;
    dip_ind = [];
    
    % This part is skipped since the dual-dipole seems to have caused an issue for Christiane. 12/06/2018 Makoto.
    %
    %     if ~isfield(STUDY.cluster(clusterIdxList(clsIdxIdx)),'dipole')
    %         STUDY = std_centroid(STUDY,ALLEEG, clusterIdxList(clsIdxIdx) , 'dipole');
    %     elseif isempty(STUDY.cluster(clusterIdxList(clsIdxIdx)).dipole)
    %         STUDY = std_centroid(STUDY,ALLEEG, clusterIdxList(clsIdxIdx) , 'dipole');
    %     end

    for k = 1:len
        abset   = STUDY.datasetinfo(STUDY.cluster(clusterIdxList(clsIdxIdx)).sets(1,k)).index;
        subject = STUDY.datasetinfo(STUDY.cluster(clusterIdxList(clsIdxIdx)).sets(1,k)).subject;
        if ~isfield(ALLEEG(abset), 'dipfit')
            warndlg2(['No dipole information available in dataset ' ALLEEG(abset).filename ' , abort plotting'], 'Aborting plot dipoles');
            return;
        end
        comp = STUDY.cluster(clusterIdxList(clsIdxIdx)).comps(k);
        cluster_dip_models(k).posxyz = ALLEEG(abset).dipfit.model(comp).posxyz;
        cluster_dip_models(k).momxyz = ALLEEG(abset).dipfit.model(comp).momxyz;
        cluster_dip_models(k).rv = ALLEEG(abset).dipfit.model(comp).rv;
        if strcmpi(ALLEEG(abset).dipfit.coordformat, 'spherical')
            if isfield(ALLEEG(abset).dipfit, 'hdmfile') %dipfit 2 spherical model
                load('-mat', ALLEEG(abset).dipfit.hdmfile);
                max_r = max(max_r, max(vol.r));
            else % old version of dipfit
                max_r = max(max_r,max(ALLEEG(abset).dipfit.vol.r));
            end
        end
        comp_to_disp{k} = [subject  ', ' 'IC' num2str(comp) ];
        if ~isempty(cluster_dip_models(k).posxyz)
            ndip = ndip +1;
            dip_ind = [dip_ind k];
        end
    end % finished going over cluster comps
    
    % Store the results.
    clusterDipoles{1,clsIdxIdx} = cluster_dip_models;
end