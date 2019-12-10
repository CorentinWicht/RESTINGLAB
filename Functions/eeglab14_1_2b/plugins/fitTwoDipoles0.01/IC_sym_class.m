% [sym_class] = IC_sym_class(EEG, 'ncomps',[1:35],'threshold', 35, 'sym_reg', 'LRR', 'plot', 'yes')
%
% Input:
% EEG = EEGLAB EEG structure
% Optional inputs:
% 'ncomps' = vector of component numbers (default = all ICs)
% 'threshold' = threshold value for "true" peak selection (default = 35)
% 'sym_reg' = symmetry region (default LRR = Large rectangular region)
%             other admitted values:
%             LSR = Large Square Region
%             SRR = Small Rectangular Region
%             SSR = Small Square Region
% 'plot' = it plots the scalp maps of the ICs identified as bilateral.
%          Admitted values: 'yes, 'no' (default = 'yes')
% Output:
% sym_class = Matrix which contains in the first line the identification 
%             number of the ICs entered in the function and in the second
%             line the results of the "classification":
%             1 = bilateral IC, 0 = unilateral IC

% Copyright (C) Catarina Piazza, Politecnico di Milano/Scientific Institute IRCCS Eugenio Medea
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

function [sym_class] = IC_sym_class(EEG, varargin)

minargs=1;  maxargs=9;

narginchk(minargs, maxargs)

arglist = finputcheck(varargin, ...
    {'ncomps'          'integer'  []  1:size(EEG.icaweights,1);...
    'threshold'        'integer'  []  35;...
    'sym_reg'          'string'   {'LRR', 'LSR', 'SRR', 'SSR'}  'LRR';...
    'plot'             'string'   {'yes','no'} 'yes';...
    });

if ischar(arglist)
    error(arglist);
end

th = arglist.threshold;
s_reg = arglist.sym_reg;
comp_indexes = arglist.ncomps;

voxel_size = [4,5];
sym_class = NaN(2,length(comp_indexes)); 
P = NaN(length(comp_indexes),4); % matrix that contains for each IC (line) the maximum and minimum values (normalized by delta_max) within each hemisphere (1st and 2nd columns, max and min of the left hemisphere)

% identification of the symmetry region
switch s_reg
    case 'LRR'
        step_c = voxel_size(1)*2;
        step_r = voxel_size(2)*4;
    case 'LSR'
        step_c = voxel_size(1)*3;
        step_r = voxel_size(2)*3;
    case 'SRR'
        step_c = voxel_size(1);
        step_r = voxel_size(2)*4;
    case 'SSR'
        step_c = voxel_size(1)*2;
        step_r = voxel_size(2)*2;
end

for n = 1:length(comp_indexes)
    comp = comp_indexes(n);
    sym_class(1,n)= comp;
    [null, map_val] = topoplot(EEG.icawinv(:, comp),EEG.chanlocs,'intrad',0.5,'numcontour',0,'noplot','on');
    
    % computation of the mean value of the map_val matrix, which is used for delta_max computation
    T_values = NaN(length(find(isnan(map_val)==0)),1); % T_values is the vector that wil be fill in with all the non NaN element of the matrix map_val
    cont = 1;
    for x = 1:size(map_val,1)
        for y = 1:size(map_val,2)
            if isnan(map_val(x,y))==0
                T_values(cont)= map_val(x,y);
                cont = cont+1;
            end
        end
    end
    m_value = mean(T_values);  % mean value of the map_val matrix
    
    T = NaN(size(map_val,1),size(map_val,2)); % matrix T will be created dividing the map_val matrix in 5x4 element voxels and assigning each voxel the mean value of its matrix elements (voxels with less than 70% of elements ~= NaN are assigned as NaN)
    
    % T_left creation
    Tleft= map_val(2:size(map_val,1)-1,2:(size(map_val,2)-1)/2); % (non vengono considerati prima colonna e prima/ultima riga perch� contengono un solo valore)
     Vleft=[]; % matrix that will be fill in with: the value of each voxel (1st column) and the position of the voxel inside the T_left matrix, line (2nd column) and column (3rd column)
    cont2 = 1;
    for i=1:voxel_size(2):(size(map_val,1)-1)-voxel_size(2)
        for j=1:voxel_size(1):((size(map_val,2)-1)/2)-voxel_size(1)
            Ttemp = Tleft(i:i+voxel_size(2)-1,j:j+voxel_size(1)-1); % contiene valori del voxel corrente
            [r,c] = find(isnan(Ttemp)==0);
            if length(r)>=(size(Ttemp,1)*size(Ttemp,2))*70/100; % (se nella roi almeno il 70% dei pixel sono valori numerici)
                v = NaN(length(r),1); % vector which contains the no NaN values of the current voxel
                for k = 1:length(v)
                    v(k)= Ttemp(r(k),c(k));
                end
                v_med = mean(v); % value assigned to the current voxel
                for k = 1:length(v)
                    Ttemp(r(k),c(k))= v_med;
                end
                Vleft(cont2,:)=[v_med,i,j];
                cont2=cont2+1;
            else
                Ttemp = NaN(voxel_size(2),voxel_size(1)); % (se ho un numero di NaN > del 30% la roi corrente assume valore NaN)
            end
            Tleft(i:i+voxel_size(2)-1,j:j+voxel_size(1)-1) = Ttemp;
        end
    end
    T(2:size(map_val,1)-1,2:(size(map_val,2)-1)/2) = Tleft;
    
    % T_right creation
    Tright= map_val(2:size(map_val,1)-1,round(size(map_val,2)/2)+1:size(map_val,2)-1);
    Vright=[];
    cont3 = 1;
    for i=1:voxel_size(2):(size(map_val,1)-1)-voxel_size(2)
        for j=1:voxel_size(1):((size(map_val,2)-1)/2)-voxel_size(1)
            Ttemp = Tright(i:i+voxel_size(2)-1,j:j+voxel_size(1)-1);
            [r,c] = find(isnan(Ttemp)==0);
            if length(r)>=(size(Ttemp,1)*size(Ttemp,2))*70/100;
                v = NaN(length(r),1);
                for k = 1:length(v)
                    v(k)= Ttemp(r(k),c(k));
                end
                v_med = mean(v);
                for k = 1:length(v)
                    Ttemp(r(k),c(k))= v_med;
                end
                Vright(cont3,:)=[v_med,i,j];
                cont3=cont3+1;
            else
                Ttemp = NaN(voxel_size(2),voxel_size(1));
            end
            Tright(i:i+voxel_size(2)-1,j:j+voxel_size(1)-1) = Ttemp;
        end
    end
    T(2:size(map_val,1)-1,round(size(map_val,2)/2)+1:size(map_val,2)-1) = Tright;
    
    % % plot of T
    % figure; surf(T);
    
    % delta_max computation
    MAX = max(max(T)); % maximum value of the whole matrix T
    MIN = min(min(T)); % minimum value of the whole matrix T
    if abs(m_value-MAX)>abs(m_value-MIN)
        delta_max = abs(m_value-MAX);
    else
        delta_max = abs(m_value-MIN);
    end
    
    % maximum and minimum values of the left hemisphere normalized by  
    % delta_max
    P(n,1)= (abs(m_value-max(Vleft(:,1)))*100)/delta_max; % p_max1;
    P(n,2)= (abs(m_value-min(Vleft(:,1)))*100)/delta_max; % p_min1;
    
    % maximum and minimum values of the right hemisphere normalized by  
    % delta_max
    P(n,3)= (abs(m_value-max(Vright(:,1)))*100)/delta_max; % p_max2;
    P(n,4)= (abs(m_value-min(Vright(:,1)))*100)/delta_max; % p_min2;
    
    picchi = find(P(n,:) >= th); %true peak identification
    % evaluation of the number and the position of true peaks 
    if isempty(picchi==1) || length(picchi)==1
        sym_class(2,n) = 0;
    else
        switch length(picchi)
            case 2
                if P(n,1) >= th && P(n,3) >= th
                    r = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),2);
                    c = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),3);
                    switch c
                        case 1
                            c = 1+voxel_size(1)*7;
                        case 1+voxel_size(1)
                            c = 1+voxel_size(1)*6;
                        case 1+voxel_size(1)*2
                            c = 1+voxel_size(1)*5;
                        case 1+voxel_size(1)*3
                            c = 1+voxel_size(1)*4;
                        case 1+voxel_size(1)*4
                            c = 1+voxel_size(1)*3;
                        case 1+voxel_size(1)*5
                            c = 1+voxel_size(1)*2;
                        case 1+voxel_size(1)*6
                            c = 1+voxel_size(1);
                        case 1+voxel_size(1)*7
                            c = 1;
                    end
                    % check of the position of the max peak in the right
                    % hemisphere (the symmetry region and the additional
                    % constraints(centrality criterion and boundary criterion)
                    % are considered)
                    r2 = Vright(find(Vright(:,1)== max(Vright(:,1))),2);
                    c2 = Vright(find(Vright(:,1)== max(Vright(:,1))),3);
                    if (r2 <= r+step_r && r2 >= r-step_r) && ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11))
                        if (c2 <= c+step_c && c2 >= c-step_c) && (c~=1 && c2~=1) &&(c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2)
                            sym_class(2,n) = 1;
                        else
                            sym_class(2,n) = 0;
                        end
                    else
                        sym_class(2,n) = 0;
                    end
                else if P(n,2) >= th && P(n,4) >= th
                        r = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),2);
                        c = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),3);
                        switch c
                            case 1
                                c = 1+voxel_size(1)*7;
                            case 1+voxel_size(1)
                                c = 1+voxel_size(1)*6;
                            case 1+voxel_size(1)*2
                                c = 1+voxel_size(1)*5;
                            case 1+voxel_size(1)*3
                                c = 1+voxel_size(1)*4;
                            case 1+voxel_size(1)*4
                                c = 1+voxel_size(1)*3;
                            case 1+voxel_size(1)*5
                                c = 1+voxel_size(1)*2;
                            case 1+voxel_size(1)*6
                                c = 1+voxel_size(1);
                            case 1+voxel_size(1)*7
                                c = 1;
                        end
                        % check of the position of the max peak in the right
                        % hemisphere (the symmetry region and the additional
                        % constraints(centrality criterion and boundary criterion)
                        % are considered)
                        r2 = Vright(find(Vright(:,1)== min(Vright(:,1))),2);
                        c2 = Vright(find(Vright(:,1)== min(Vright(:,1))),3);
                        if (r2 <= r+step_r && r2 >= r-step_r) && ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11))
                            if (c2 <= c+step_c && c2 >= c-step_c) && (c~=1 && c2~=1) &&(c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2)
                                
                                sym_class(2,n) = 1;
                            else
                                sym_class(2,n) = 0;
                            end
                        else
                            sym_class(2,n) = 0;
                        end
                    else
                        sym_class(2,n)=0;
                    end
                end
            case 3
                if P(n,1) >= th && P(n,3) >= th
                    r = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),2);
                    c = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),3);
                    switch c
                        case 1
                            c = 1+voxel_size(1)*7;
                        case 1+voxel_size(1)
                            c = 1+voxel_size(1)*6;
                        case 1+voxel_size(1)*2
                            c = 1+voxel_size(1)*5;
                        case 1+voxel_size(1)*3
                            c = 1+voxel_size(1)*4;
                        case 1+voxel_size(1)*4
                            c = 1+voxel_size(1)*3;
                        case 1+voxel_size(1)*5
                            c = 1+voxel_size(1)*2;
                        case 1+voxel_size(1)*6
                            c = 1+voxel_size(1);
                        case 1+voxel_size(1)*7
                            c = 1;
                    end
                    % check of the position of the max peak in the right
                    % hemisphere (the symmetry region and the additional
                    % constraints(centrality criterion and boundary criterion)
                    % are considered)
                    r2 = Vright(find(Vright(:,1)== max(Vright(:,1))),2);
                    c2 = Vright(find(Vright(:,1)== max(Vright(:,1))),3);
                    if (r2 <= r+step_r && r2 >= r-step_r) && ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11))
                        if (c2 <= c+step_c && c2 >= c-step_c) && (c~=1 && c2~=1) &&(c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2)
                            sym_class(2,n) = 1;
                        else
                            sym_class(2,n) = 0;
                        end
                    else
                        sym_class(2,n) = 0;
                    end
                else if P(n,2) >= th && P(n,4) >= th
                        r = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),2);
                        c = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),3);
                        switch c
                            case 1
                                c = 1+voxel_size(1)*7;
                            case 1+voxel_size(1)
                                c = 1+voxel_size(1)*6;
                            case 1+voxel_size(1)*2
                                c = 1+voxel_size(1)*5;
                            case 1+voxel_size(1)*3
                                c = 1+voxel_size(1)*4;
                            case 1+voxel_size(1)*4
                                c = 1+voxel_size(1)*3;
                            case 1+voxel_size(1)*5
                                c = 1+voxel_size(1)*2;
                            case 1+voxel_size(1)*6
                                c = 1+voxel_size(1);
                            case 1+voxel_size(1)*7
                                c = 1;
                        end
                        % check of the position of the max peak in the right
                        % hemisphere (the symmetry region and the additional
                        % constraints(centrality criterion and boundary criterion)
                        % are considered)
                        r2 = Vright(find(Vright(:,1)== min(Vright(:,1))),2);
                        c2 = Vright(find(Vright(:,1)== min(Vright(:,1))),3);
                        if (r2 <= r+step_r && r2 >= r-step_r) && ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11))
                            if (c2 <= c+step_c && c2 >= c-step_c) && (c~=1 && c2~=1) &&(c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2)
                                sym_class(2,n) = 1;
                            else
                                sym_class(2,n) = 0;
                            end
                        else
                            sym_class(2,n) = 0;
                        end
                    else
                        sym_class(2,n)=0;
                    end
                end
            case 4
                who = find(P(n,:)== max(P(n,:)));
                if who == 1 || who == 3
                    r = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),2);
                    c = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),3);
                    switch c
                        case 1
                            c = 1+voxel_size(1)*7;
                        case 1+voxel_size(1)
                            c = 1+voxel_size(1)*6;
                        case 1+voxel_size(1)*2
                            c = 1+voxel_size(1)*5;
                        case 1+voxel_size(1)*3
                            c = 1+voxel_size(1)*4;
                        case 1+voxel_size(1)*4
                            c = 1+voxel_size(1)*3;
                        case 1+voxel_size(1)*5
                            c = 1+voxel_size(1)*2;
                        case 1+voxel_size(1)*6
                            c = 1+voxel_size(1);
                        case 1+voxel_size(1)*7
                            c = 1;
                    end
                    % check of the position of the max peak in the right
                    % hemisphere (the symmetry region and the additional
                    % constraints(centrality criterion and boundary criterion)
                    % are considered)
                    r2 = Vright(find(Vright(:,1)== max(Vright(:,1))),2);
                    c2 = Vright(find(Vright(:,1)== max(Vright(:,1))),3);
                    if ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11)) && (c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2) && (c~=1 && c2~=1)
                        if (r2 <= r+step_r && r2 >= r-step_r)
                            if (c2 <= c+step_c && c2 >= c-step_c)
                                sym_class(2,n) = 1;
                            else
                                sym_class(2,n) = 0;
                            end
                        else
                            sym_class(2,n) = 0;
                        end
                    else
                        r = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),2);
                        c = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),3);
                        switch c
                            case 1
                                c = 1+voxel_size(1)*7;
                            case 1+voxel_size(1)
                                c = 1+voxel_size(1)*6;
                            case 1+voxel_size(1)*2
                                c = 1+voxel_size(1)*5;
                            case 1+voxel_size(1)*3
                                c = 1+voxel_size(1)*4;
                            case 1+voxel_size(1)*4
                                c = 1+voxel_size(1)*3;
                            case 1+voxel_size(1)*5
                                c = 1+voxel_size(1)*2;
                            case 1+voxel_size(1)*6
                                c = 1+voxel_size(1);
                            case 1+voxel_size(1)*7
                                c = 1;
                        end
                        % check of the position of the max peak in the right
                        % hemisphere (the symmetry region and the additional
                        % constraints(centrality criterion and boundary criterion)
                        % are considered)
                        r2 = Vright(find(Vright(:,1)== min(Vright(:,1))),2);
                        c2 = Vright(find(Vright(:,1)== min(Vright(:,1))),3);
                        if (r2 <= r+step_r && r2 >= r-step_r) && ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11))
                            if (c2 <= c+step_c && c2 >= c-step_c) && (c~=1 && c2~=1) &&(c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2)
                                sym_class(2,n) = 1;
                            else
                                sym_class(2,n) = 0;
                            end
                        else
                            sym_class(2,n) = 0;
                        end
                    end
                else
                    r = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),2);
                    c = Vleft(find(Vleft(:,1)== min(Vleft(:,1))),3);
                    switch c
                        case 1
                            c = 1+voxel_size(1)*7;
                        case 1+voxel_size(1)
                            c = 1+voxel_size(1)*6;
                        case 1+voxel_size(1)*2
                            c = 1+voxel_size(1)*5;
                        case 1+voxel_size(1)*3
                            c = 1+voxel_size(1)*4;
                        case 1+voxel_size(1)*4
                            c = 1+voxel_size(1)*3;
                        case 1+voxel_size(1)*5
                            c = 1+voxel_size(1)*2;
                        case 1+voxel_size(1)*6
                            c = 1+voxel_size(1);
                        case 1+voxel_size(1)*7
                            c = 1;
                    end
                    % check of the position of the max peak in the right
                    % hemisphere (the symmetry region and the additional
                    % constraints(centrality criterion and boundary criterion)
                    % are considered)
                    r2 = Vright(find(Vright(:,1)== min(Vright(:,1))),2);
                    c2 = Vright(find(Vright(:,1)== min(Vright(:,1))),3);
                    if ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11)) && (c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2) && (c~=1 && c2~=1)
                        if (r2 <= r+step_r && r2 >= r-step_r)
                            if (c2 <= c+step_c && c2 >= c-step_c)
                                sym_class(2,n) = 1;
                            else
                                sym_class(2,n) = 0;
                            end
                        else
                            sym_class(2,n) = 0;
                        end
                    else
                        r = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),2);
                        c = Vleft(find(Vleft(:,1)== max(Vleft(:,1))),3);
                        switch c
                            case 1
                                c = 1+voxel_size(1)*7;
                            case 1+voxel_size(1)
                                c = 1+voxel_size(1)*6;
                            case 1+voxel_size(1)*2
                                c = 1+voxel_size(1)*5;
                            case 1+voxel_size(1)*3
                                c = 1+voxel_size(1)*4;
                            case 1+voxel_size(1)*4
                                c = 1+voxel_size(1)*3;
                            case 1+voxel_size(1)*5
                                c = 1+voxel_size(1)*2;
                            case 1+voxel_size(1)*6
                                c = 1+voxel_size(1);
                            case 1+voxel_size(1)*7
                                c = 1;
                        end
                        % check of the position of the max peak in the right
                        % hemisphere (the symmetry region and the additional
                        % constraints(centrality criterion and boundary criterion)
                        % are considered)
                        r2 = Vright(find(Vright(:,1)== max(Vright(:,1))),2);
                        c2 = Vright(find(Vright(:,1)== max(Vright(:,1))),3);
                        if (r2 <= r+step_r && r2 >= r-step_r) && ((voxel_size(2)*2<r && r<=voxel_size(2)*11) || (voxel_size(2)*2<r2 && r2<=voxel_size(2)*11))
                            if (c2 <= c+step_c && c2 >= c-step_c) && (c~=1 && c2~=1) &&(c>=1+voxel_size(1)*2 || c2>=1+voxel_size(1)*2)
                                sym_class(2,n) = 1;
                            else
                                sym_class(2,n) = 0;
                            end
                        else
                            sym_class(2,n) = 0;
                        end
                    end
                end
        end
    end
end

% plot of bilateral ICs
if strcmp(arglist.plot,'yes')==1
    plot_comp = [];
    plot_comp_index = find(sym_class(2,:)==1);
    if isempty(plot_comp_index)== 0  
        for index = 1:length(plot_comp_index)
            plot_comp(index) = sym_class(1,plot_comp_index(index));
        end
        pop_selectcomps(EEG, plot_comp);
    else
        fprintf('\nNo bilateral components identified\n');
    end
end