function MI_EDS_combined = climada_measures_impact_report(measures_impact,report_xls_file)
% climada_measures_impact_report
% MODULE:
%   climada core
% NAME:
%   climada_measures_impact_report
% PURPOSE:
%   combine EDS's of measures_impact structure array, then print to excel
% CALLING SEQUENCE:
%   MI_EDS_combined=climada_measures_impact_report(measures_impact,report_xls_file)
% EXAMPLE:
%   shape_plotter(shapes,'attribute name','linewidth',2,'color','r')
% INPUTS: 
%   measures_impact:   measures impact struct, see climada_measures_impact
% OPTIONAL INPUT PARAMETERS:
%   report_xls_file:   xls file to save to
% OUTPUTS:
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150623, init
%-

if ~exist('measures_impact'     ,'var'),    measures_impact = [];   end
if ~exist('report_xls_file'     ,'var'),    report_xls_file = '';   end
if isempty(measures_impact)            ,    return;                 end

sheet = num2str(measures_impact(1).EDS(1).reference_year);
MI_EDS_combined = measures_impact(1).EDS; % control init

fprintf('measures impact report for:\n')
for k=1:length(MI_EDS_combined)
    MI_EDS_combined(k).ED_at_centroid   = zeros(size(MI_EDS_combined(k).ED_at_centroid)); % init
    fld(k).MI_at_centroid               = zeros(size(MI_EDS_combined(k).ED_at_centroid)); % init
    MI_EDS_combined(k).peril_ID         = '';
    if ~strcmp(MI_EDS_combined(k).annotation_name,'control')
        fprintf('\t%s\n',MI_EDS_combined(k).annotation_name)
    end
end

for i = 1:length(measures_impact) % loop through each hazard/entity combo
    % ED at centroid w/o measure
    ctrl_ED_at_c = measures_impact(i).EDS(end).ED_at_centroid;
    
    missing = ~ismember({measures_impact(i).EDS(:).annotation_name},{MI_EDS_combined(:).annotation_name});
    if any(missing)
        for m = 1:length(measures_impact(i).EDS)
            if missing(m)
                fprintf('\t%s\n',measures_impact(i).EDS(m).annotation_name)
            end
        end
        % add new measure to struct
        MI_EDS_combined             = [MI_EDS_combined measures_impact(i).EDS(missing)];
        fld(end+1).MI_at_centroid   = zeros(size(fld(1).MI_at_centroid)); % init for new measure
    end
    
    for j = 1:length(MI_EDS_combined)
        % for each measure already in MI_EDS_combined, see if it exists in
        % the next measures_impact struct, add to MI_EDS_combined if it
        % does, else add entry to MI_EDS_combined
        ndx = find(strcmpi(MI_EDS_combined(j).annotation_name,{measures_impact(i).EDS(:).annotation_name}));
        if ~isempty(ndx)
            MI_EDS_combined(j).ED_at_centroid  	= MI_EDS_combined(j).ED_at_centroid + measures_impact(i).EDS(ndx).ED_at_centroid;
            fld(j).MI_at_centroid               = fld(j).MI_at_centroid + (ctrl_ED_at_c -measures_impact(i).EDS(ndx).ED_at_centroid);
            
            MI_EDS_combined(j).hazard(end+1)    = measures_impact(i).EDS(ndx).hazard;
            MI_EDS_combined(j).peril_ID         = [MI_EDS_combined(j).peril_ID measures_impact(i).EDS(ndx).peril_ID ' | '];
        else
            MI_EDS_combined(j).hazard           = measures_impact(i).EDS(ndx).hazard;
            MI_EDS_combined(j).peril_ID         = measures_impact(i).EDS(ndx).peril_ID;

            MI_EDS_combined(j).ED_at_centroid   = MI_EDS_combined(j).ED_at_centroid+ctrl_ED_at_c;
            fld(j).MI_at_centroid               = fld(j).MI_at_centroid + 0;
        end
    end
end
% assign fld.MI_at_centroid to MI_EDS_combined
for l = 1:length(MI_EDS_combined)
    MI_EDS_combined(l).MI_at_centroid = fld(l).MI_at_centroid;
end

% move control to end if it isn't already
for n = 1:length(MI_EDS_combined)
    if strcmpi(MI_EDS_combined(n).annotation_name,'control')
        control_ndx = n;
        break
    end
end
MI_EDS_combined(end+1) = MI_EDS_combined(n); % add control to end
MI_EDS_combined(n)     = [];                 % delete original control
    
% save(or not)
if ~isempty(report_xls_file)
    if ~strcmp(report_xls_file,'NO_SAVE')
        climada_EDS_ED_at_centroid_report_xls(MI_EDS_combined,EDS_report_xls,sheet,'MI_at_centroid');
    end
else
    report_xls_file=[module_data_dir filesep 'entities' filesep '*.mat'];
    [fN, pN] = uiputfile(report_xls_file, 'Save measures imnpact report as:');
    if isequal(fN,0) || isequal(pN,0)
        return; % cancel
    else
        climada_EDS_ED_at_centroid_report_xls(MI_EDS_combined,EDS_report_xls,sheet,'MI_at_centroid');
    end
end
