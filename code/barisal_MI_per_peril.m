function barisal_MI_per_peril(measures_impact5,measures,peril_IDs)
% climada barisal demo
% NAME:
%   barisal_EDS_DFC, special version for BARISAL
% PURPOSE:
%   create xls-report for measures impact per peril (and overall effectiveness 
%   of measures)
%
%-

fig=[];legend_str=[];

global climada_global
if ~climada_init_vars,return;end % init/import global variables

measure_names = unique(measures.name);
no_measures = numel(unique(measures.name))+1;
matr = cell(no_measures+2,length(peril_IDs)+1);
matr{1,1} = 'Scenario 1, 2050, moderate climate change';
matr{1,2} = 'Annual expected damage (AED, BDT) per peril type, without measures (baseline) and with measures';
matr{2+1,1} = 'Baseline';
matr(2+2:3+no_measures-1,1) = measure_names;
for p_i = 5:8%1:length(peril_IDs)
    
    %peril type
    if numel(measures_impact5(p_i).EDS(1).event_ID) == 4
        peril_name = 'cyclone';
    elseif numel(measures_impact5(p_i).EDS(1).event_ID) == 29
        peril_name = 'monsoon';
    else
        peril_name = 'cyclone';
    end
        
    switch measures_impact5(p_i).EDS(1).peril_ID
        case 'FL'
            peril_name = sprintf('Remaining AED due to flood depth %s', peril_name);
        case 'FL_'
            peril_name = sprintf('Remaining AED due to flood duration %s', peril_name);
        case 'TC'
            peril_name = 'Remaining AED due to cyclone wind';
    end
    matr{2,p_i+3} = peril_name;
    
    % baseline AED
    matr{2+1,p_i+3} = measures_impact5(p_i).ED(end);
    
    for m_i = 1:no_measures-1
        %find measure in measures_impact.EDS structure
        indx = find(strcmp(measure_names(m_i),{measures_impact5(p_i).EDS.annotation_name}));
        if ~isempty(indx)
            indx = indx(1);
            matr{3+m_i,p_i+3} = measures_impact5(p_i).ED(indx);
        else
            % baseline value
            matr{3+m_i,p_i+3} = measures_impact5(p_i).ED(end);
        end
    end
end
%total AED
matr{2,2} = 'Total AED (BDT)';
matr{2,3} = 'Overall effectiveness (averted damage in %)';
for m_i = 1:no_measures
    matr{2+m_i,2} =sum([matr{2+m_i,3:end}]); 
    % benefit in percentage
    matr{2+m_i,3} = (matr{2+1,2}-matr{2+m_i,2})/matr{2+1,2}; 
end

barisal_data_dir= [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
results_dir = [barisal_data_dir filesep 'results'];
xls_file = [results_dir filesep 'Measures_impact_per_peril_type_Scenario1_2050.xls'];
xlswrite(xls_file,matr)




