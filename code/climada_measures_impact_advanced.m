function measures_impact=climada_measures_impact_advanced(entity,hazard,measures_impact_reference,measures,map_risk_premium,sanity_check)
% special measures impact function for Barisal
% MODULE:
%   barisal_demo
% NAME:
%   climada_measures_impact_advanced
% PURPOSE:
%   special measures impact function for barisal, that includes hazard
%   modification files (i.e. absolute or percentage reduction of flood
%   depths and duration)
%   see climada_measures_impact for more information
% CALLING SEQUENCE:
%   measures_impact = climada_measures_impact_advanced(entity,hazard,measures_impact_reference,measures,map_risk_premium,sanity_check)
% EXAMPLE:
%   measures_impact = climada_measures_impact_advanced % all prompted for
%   hazard_set_file='...\climada\data\hazards\TCNA_A_Probabilistic.mat';
%   measures_impact=climada_measures_impact(climada_entity_read('',hazard_set_file),hazard_set_file,'no')
%   measures_impact=climada_measures_impact('','','','',1) % all interactive, show risk premium map
% INPUTS:
%   entity: a read and encoded assets and damagefunctions file, see climada_assets_encode(climada_assets_read)
%       > promted for if not given
%   hazard: either a hazard set (struct) or a hazard set file (.mat with a struct)
%       > promted for if not given
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   measures_impact: a structure with
%       EDS(measure_i): the event damage set for each measure, last one EDS(end) for no measures
%       ED(measure_i): the annual expected damage to the assets under measure_i,
%           last one ED(end) for no measures
%       benefit(measure_i): the benefit of measure_i
%       cb_ratio(measure_i): the cost/benefit ratio of measure_i
%       measures: just a copy of measures, so we have all we need together
%       title_str: a meaningful title, of the format: measures @ assets | hazard
%       NOTE: currently measures_impact is also stored (with a lengthy filename)
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150611, created from
%               original climada_measures_impact; added functionality for
%               hazard_event_set_operator & measures_distributed
% Gilles Stassen, gillesstassen@hotmail.com, 20150616, entity switch added
% Lea Mueller, muellele@gmail.com, 20150902, rename to hazard_intensity_impact_b from hazard_intensity_impact
% Lea Mueller, muellele@gmail.com, 20151116, add regional scope, add documentation
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

measures_impact=[];

% poor man's version to check arguments
if ~exist('entity','var'),entity=[];end
if ~exist('hazard','var'),hazard=[];end
if ~exist('measures_impact_reference','var'),measures_impact_reference=[];end
if ~exist('measures','var'),measures='';end
if ~exist('map_risk_premium','var'),map_risk_premium=0;end
if ~exist('sanity_check','var'),sanity_check=1;end

% PARAMETERS

% prompt for entity if not given
if isempty(entity) % local GUI
    entity=[climada_global.data_dir filesep 'entities' filesep '*.mat'];
    [filename, pathname] = uigetfile(entity, 'Select encoded entity:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        entity=fullfile(pathname,filename);
    end
end
% load the entity, if a filename has been passed
if ~isstruct(entity)
    entity_file=entity;entity=[];
    vars = whos('-file', entity_file);
    load(entity_file);
    if ~strcmp(vars.name,'entity')
        entity = eval(vars.name);
        clear (vars.name)
    end
end

% prompt for hazard if not given
if isempty(hazard) % local GUI
    hazard=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard event set for EDS calculation:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard=fullfile(pathname,filename);
    end
end
% load the hazard, if a filename has been passed
if ~isstruct(hazard)
    hazard_file=hazard;hazard=[];
    vars = whos('-file', hazard_file);
    load(hazard_file);
    if ~strcmp(vars.name,'hazard')
        hazard = eval(vars.name);
        clear (vars.name)
    end
end
hazard=climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files

% prompt for reference result if not given
if isempty(measures_impact_reference)
    measures_impact_reference=[climada_global.data_dir filesep 'results' filesep '*.mat'];
    [filename, pathname] = uigetfile(measures_impact_reference, 'Select reference results (skip if not future):');
    if isequal(filename,0) || isequal(pathname,0)
        measures_impact_reference=''; % re-set to empty
    else
        measures_impact_reference = fullfile(pathname,filename);
    end
end

if ~isempty(measures_impact_reference)
    % load if a filename has been passed
    if ~isstruct(measures_impact_reference)
        if strcmp(measures_impact_reference,'no')
            measures_impact_reference = '';
        else
            measures_impact_reference_file = measures_impact_reference; measures_impact_reference = [];
            vars = whos('-file', measures_impact_reference_file);
            load(measures_impact_reference_file);
            if ~strcmp(vars.name,'measures_impact_reference')
                measures_impact_reference = eval(vars.name);
                clear (vars.name)
            end
        end
    end
    if ~isempty(measures_impact_reference)
        % reference is aleays today, warn if not
        reference_year = measures_impact_reference.EDS(end).reference_year;
        if reference_year ~= climada_global.present_reference_year
            %fprintf('WARNING: reference year for reference results is %i (should be %i)\n',reference_year,climada_global.present_reference_year);
        end
    end
end

if isfield(entity,'measures')
    if isempty(measures)
        measures = entity.measures;
    end
elseif isempty(measures)
    measures = 'ASK';
end
% prompt for measures if not given
if ~isstruct(measures)
    if strcmp(measures,'ASK')
        measures = [climada_global.data_dir filesep 'measures' filesep '*.mat'];
        [filename, pathname] = uigetfile(measures, 'Select measures:');
        if isequal(filename,0) || isequal(pathname,0)
            return; % cancel
        else
            measures = fullfile(pathname,filename);
        end
    end
end
% load the measures, if a filename has been passed
if ~isstruct(measures)
    measures_file=measures;measures=[];
    vars = whos('-file', measures_file);
    load(measures_file);
    if ~strcmp(vars.name,'measures')
        measures = eval(vars.name);
        clear (vars.name)
    end
end

% encode measures
measures = climada_measures_encode(measures);

% reduce measures struct to contain only those with the right peril ID
if isfield(measures,'peril_ID') && isfield(hazard,'peril_ID')
    measures_ori = measures; % backup copy
    m_flds = fieldnames(measures);
    rm_ndx = [];
    for measure_i = 1:length(measures.cost)
        if ~strcmp(measures.peril_ID{measure_i},hazard.peril_ID) && ~isempty(measures.peril_ID{measure_i})
            rm_ndx = [rm_ndx 1]; %peril_ID does not match, therefore remove this measure
        else
            rm_ndx = [rm_ndx 0]; % keep this measures
        end
        
        % some hardwiring for Barisal
        % 20151117, add flood resilient buildings (project package)
        if ismember(lower(measures.name{measure_i}),...
                   {'khals deepened 1m' 'khals widened 4m' 'ponds retain' ...
                    'drainage 5cm less water logging' 'drainage 15cm less water logging' ...
                    'solid waste management' 'pervious pavement'}) ... %'FL depth project package' 'flood resilient buildings' 'FL duration project package'}) 
                    ... 
                && ~isempty(strfind(hazard.filename,'cyclone'))
            rm_ndx(end) = 1;
        end
    end
    measures.color_RGB = measures.color_RGB(~rm_ndx,:);
    for m_fld_i = 1:length(m_flds)
        if length(measures.(m_flds{m_fld_i})) == length(measures_ori.cost)
            measures.(m_flds{m_fld_i}) = measures.(m_flds{m_fld_i})(~rm_ndx);
        end
    end
end
if all(rm_ndx)
    cprintf([1 0.5 0],'WARNING: no measures found with peril ID %s\n',hazard.peril_ID);
    measures = climada_measures_construct([],1); % neutral measure
    measures.name{1} = 'control';
%     return;
end

% check for correct encoding to the hazard
if ~isfield(entity.assets,'hazard')
    force_re_encode=1; % entity not encoded yet
else
    % assets have been encoded, check whether for same hazard
    if ~strcmp(entity.assets.hazard.filename,hazard.filename)
        % not the same hazard set, force re-encoding to be on the safe side
        % (it's likely that the climate scenario hazard event set is valid at
        % the exact same centroids, but encoding is much faster than any later
        % troubleshooting ;-)
        force_re_encode=1;
    else
        force_re_encode=0; % entity encoded
    end
end

if force_re_encode
    entity=climada_assets_encode(entity,hazard);
end

% loop over all measures and calculate the corresponding EDSs
% -----------------------------------------------------------

% make (backup) copies of the original data in entity:
entity_orig_damagefunctions = entity.damagefunctions;
orig_assets_DamageFunID     = entity.assets.DamageFunID;
orig_assets = entity.assets; % backup

if isfield(measures,'damagefunctions') % append measure's vulnerabilities
    entity.damagefunctions.DamageFunID = measures.damagefunctions.DamageFunID;
    entity.damagefunctions.Intensity   = measures.damagefunctions.Intensity;
    entity.damagefunctions.MDD         = measures.damagefunctions.MDD;
    entity.damagefunctions.PAA         = measures.damagefunctions.PAA;
    
    %entity.damagefunctions.DamageFunID = [entity.damagefunctions.DamageFunID,measures.damagefunctions.DamageFunID];
    %entity.damagefunctions.Intensity   = [entity.damagefunctions.Intensity,measures.damagefunctions.Intensity];
    %entity.damagefunctions.MDD         = [entity.damagefunctions.MDD,measures.damagefunctions.MDD];
    %entity.damagefunctions.PAA         = [entity.damagefunctions.PAA,measures.damagefunctions.PAA];
end

orig_damagefunctions = entity.damagefunctions; % copy of this table for speedup
n_measures           = length(measures.cost);

% check measures, correct regional_scope if necessary to ensure size is aligned with number of assets
measures = climada_measures_check(measures,entity.assets);

%fprintf('assessing impacts of %i measures:\n',n_measures);

ED_risk_transfer = zeros(1,n_measures+1); % allocate
hazard_switched = 0; % indicated whether a special hazard set for a given measure is used
entity_switched = 0; % indicated whether a special entity for a given measure is used
assets_switched = 0; % indicated whether a special assets for a given measure is used

% calculate control (with no measure) right from the beginning
EDS(n_measures+1) = climada_EDS_calc(entity,hazard,'control','',0,sanity_check);

for measure_i = 1:n_measures %+1 % last with no measures

    %if measure_i <= n_measures
        
    if strcmpi(measures.name{measure_i},'Flood resilient crops')
        pause(0);
    end

    %fprintf('assessing impact of measure %i\n',measure_i); % TEST

    % special treatment if an alternate hazard set is given
    if isfield(measures,'hazard_event_set')
        measures_hazard_set_name = measures.hazard_event_set{measure_i};
        if ~strcmp(measures_hazard_set_name,'nil')
            orig_hazard = hazard;
            if ~exist(measures_hazard_set_name,'file')
                % only filename given in measures tab, add path:
                if exist('hazard_file','var')
                    [hazard_dir] = fileparts(hazard_file);
                else
                    hazard_dir = [climada_global.data_dir filesep 'hazards']; % default
                end
                measures_hazard_file = [hazard_dir filesep measures_hazard_set_name];
            else
                % contains path already
                measures_hazard_file = measures_hazard_set_name;
            end
            [fP,fN,fE] = fileparts(measures_hazard_file);
            if isempty(fE),measures_hazard_file = [fP filesep fN '.mat'];end % append .mat
            if exist(measures_hazard_file,'file')
                % orig_hazard = hazard; % backup
                cprintf([0 0 1],'NOTE: measure %i, modified hazard according to %s\n',measure_i,measures_hazard_file);
                %load(measures_hazard_file);                    

                if isfield(measures,'hazard_event_set_operator')
                    switch char(measures.hazard_event_set_operator{measure_i}(1))
                        case 'p',  	op = @plus;
                        case 'm', 	op = @minus;
                        case 't', 	op = @times;
                        case 'd',   op = @divide;
                        otherwise
                            cprintf([1 0.5 0],'WARNING: operator ''%s'' not recognised\n',measures.hazard_event_set_operator{measure_i})
                            op = []; % distributed_measures will try to figure it out from the data
                    end
                    if ~strcmp(measures.hazard_event_set_operator,'nil')
                        hazard = climada_distributed_measures(measures_hazard_file,orig_hazard,op,0);
                    else
                        load(measures_hazard_file);
                    end
                else
                    load(measures_hazard_file);
                end

                hazard_switched = 1;
                % some basic consistency checks to make sure switched
                % hazard is reasonable
                try
                    hazard_checksum = abs(sum(hazard.centroid_ID-orig_hazard.centroid_ID + ...
                        hazard.lon-orig_hazard.lon + ...
                        hazard.lat-orig_hazard.lat));
                    if hazard_checksum>0
                        fprintf('WARNING: hazard might not be fully compatible with encoded assets\n');
                    end
                catch
                    % if array dimension do not match, above check fails
                    fprintf('WARNING: hazard might not at all be compatible with encoded assets\n');
                end % try
                % hazard=climada_hazard2octave(hazard); % Octave compatibility for -v7.3 mat-files
                % if abs(numel(hazard.intensity)-numel(hazard.intensity))
                %     fprintf('WARNING: hazard might not be fully compatible with EDS\n');
                % end
            else
                cprintf([1 0 0],'ERROR: measure %i, hazard NOT switched, hazard set %s not found\n',measure_i,measures_hazard_file);
            end
        else
            hazard_switched = 0; % no hazard switched
        end % measures_hazard_set_name
    end % isfield(measures,'hazard_event_set')

    % special treatment if an alternate entity is given
    if isfield(measures,'entity_file')
        measures_entity_name = measures.entity_file{measure_i};
        if ~strcmp(measures_entity_name,'nil')
            orig_entity = entity;
            if ~exist(measures_entity_name,'file')
                % only filename given in measures tab, add path:
                if exist('entity_file','var')
                    [entity_dir] = fileparts(entity_file);
                else
                    entity_dir = [climada_global.data_dir filesep 'entities']; % default
                end
                measures_entity_file = [entity_dir filesep measures_entity_name];
            else
                % contains path already
                measures_entity_file = measures_entity_name;
            end
            [fP,fN,fE] = fileparts(measures_entity_file);
            if isempty(fE),measures_entity_file = [fP filesep fN '.mat'];end % append .mat
            if exist(measures_entity_file,'file')
                cprintf([0 0 1],'NOTE: measure %i, switched entity according to %s\n',measure_i,measures_entity_file);
                load(measures_entity_file);
                entity = climada_assets_encode(entity,hazard);
                entity_switched = 1;
            else
                cprintf([1 0 0],'ERROR: measure %i, entity NOT switched, entity %s not found\n',measure_i,measures_entity_file);
            end
        else
            entity_switched = 0; % no entity switched
        end % measures_entity_name
    end % isfield(measures,'entity_file')
        
    % switch assets if needed (defined in measures.assets_file)
    assets = climada_measures_assets_change(measures,measure_i);
    if ~isempty(assets)
        entity.assets = climada_assets_encode(assets,hazard);
        assets_switched = 1;
        %save entity.assets so that it can be used in climada_EDS_ED_report
        if ~isempty(strfind(entity.assets.filename,'.mat')),save(entity.assets.filename,'entity'),end
    end

    for map_i = 1:length(measures.damagefunctions_mapping(measure_i).map_from)
        % damagefunctions mapping
        pos = orig_assets_DamageFunID==measures.damagefunctions_mapping(measure_i).map_from(map_i);
        entity.assets.DamageFunID(pos) = measures.damagefunctions_mapping(measure_i).map_to(map_i);
        %fprintf('mapping DamageFunID %i to %i ',...
        %    measures.damagefunctions_mapping(measure_i).map_from(map_i),...
        %    measures.damagefunctions_mapping(measure_i).map_to(map_i));
    end % map_i

    %if measure_i==debug_measure_i,plot(orig_damagefunctions.Intensity,orig_damagefunctions.MDD.*orig_damagefunctions.PAA,'--b');end

    entity.damagefunctions.Intensity = max(orig_damagefunctions.Intensity - measures.hazard_intensity_impact_b(measure_i),0);
    entity.damagefunctions.MDD       = max(orig_damagefunctions.MDD*measures.MDD_impact_a(measure_i)+measures.MDD_impact_b(measure_i),0);
    entity.damagefunctions.PAA       = max(orig_damagefunctions.PAA*measures.PAA_impact_a(measure_i)+measures.PAA_impact_b(measure_i),0);
    annotation_name                  = measures.name{measure_i};

    if isfield(measures,'hazard_intensity_impact_a')
        % for BACKWARD COMPATIBILITY: the following line was used until 20150907
        %%entity.damagefunctions.Intensity = max(entity.damagefunctions.Intensity*(1-measures.hazard_intensity_impact_a(measure_i)),0);
        if ~(measures.hazard_intensity_impact_a(measure_i)==0)
            entity.damagefunctions.Intensity = max(entity.damagefunctions.Intensity/measures.hazard_intensity_impact_a(measure_i),0);
        end
    end

    % add regional scope of a measure, correct .ED_at_centroid, .ED . Do not use .damage afterwards, as this holds an incorrect number
    if isfield(measures,'regional_scope')
        within_scope = measures.regional_scope(:,measure_i);
        if any(~within_scope) %measure does not affect all assets
            entity.assets.Value(~within_scope) = 0;
            assets_switched = 1;
            fprintf('Measure %d does not affect all assets\n', measure_i)
        end
    end % regional_scope
   
    % finally do the damage calculations
    EDS(measure_i) = climada_EDS_calc(entity,hazard,annotation_name,'',0,sanity_check);


%     % special routine for resilient buildings barisal in measures package
%     if strcmp(hazard.peril_ID,'FL') && isfield(entity.assets,'most_vuln_building') &&...
%             entity.assets.reference_year > climada_global.present_reference_year 
%         ent_tmp     = entity; % backup
%         test_mode   = 1;
%         switch entity.assets.reference_year
%             case 2030,  mod_A = 0.5;    mod_B = 0.3;
%             case 2050,  mod_A = 1.0;    mod_B = 0.3;
%         end
%         
%         ndx_A = strcmp('A',ent_tmp.assets.most_vuln_building);
%         ndx_B = strcmp('B',ent_tmp.assets.most_vuln_building);
%         
%         damfunIDs_A = ent_tmp.assets.DamageFunID(ndx_A);
%         damfunIDs_B = ent_tmp.assets.DamageFunID(ndx_B);
%         
%         damfun_mod_A = ismember(ent_tmp.damagefunctions.DamageFunID,damfunIDs_A) & (ent_tmp.damagefunctions.Intensity ~= 0); % don't shift the zero
%         damfun_mod_B = ismember(ent_tmp.damagefunctions.DamageFunID,damfunIDs_B) & (ent_tmp.damagefunctions.Intensity ~= 0);
%         
%         ent_tmp.damagefunctions.Intensity(damfun_mod_A) = ent_tmp.damagefunctions.Intensity(damfun_mod_A)+mod_A;
%         ent_tmp.damagefunctions.Intensity(damfun_mod_B) = ent_tmp.damagefunctions.Intensity(damfun_mod_B)+mod_B;
%         
%         if test_mode
%             figure; hold on
%             plot(entity.damagefunctions.Intensity(damfun_mod_A),entity.damagefunctions.MDD(damfun_mod_A),'r')
%             plot(entity.damagefunctions.Intensity(damfun_mod_B),entity.damagefunctions.MDD(damfun_mod_B),'b')
%             plot(ent_tmp.damagefunctions.Intensity(damfun_mod_A),ent_tmp.damagefunctions.MDD(damfun_mod_A),'g')
%             plot(ent_tmp.damagefunctions.Intensity(damfun_mod_B),ent_tmp.damagefunctions.MDD(damfun_mod_B),'m')
%             hold off
%         end
%         EDS_tmp = climada_EDS_calc(ent_tmp,hazard,annotation_name);
%         EDS(measure_i).ED_at_centroid(ndx_A | ndx_B) = EDS_tmp.ED_at_centroid(ndx_A | ndx_B);
%         EDS(measure_i).ED = sum(EDS(measure_i).ED_at_centroid);
%         EDS(measure_i).comment = 'modified for flood resilient buildings, .damage field may contain incorrect values';
%         
%         entity.damagefunctions      = orig_damagefunctions;
%     end
        
    entity.assets.DamageFunID   = orig_assets_DamageFunID; % reset damagefunctions mapping

        
    % add regional scope of a measure, correct .ED_at_centroid, .ED . Do not use .damage afterwards, as this holds an incorrect number
    if isfield(measures,'regional_scope')
        within_scope = measures.regional_scope(:,measure_i);
        if any(~within_scope) %measure does not affect all assets
            % SPECIAL CASE FOR BARISAL: 
            % overwrite assets that are not affected by this measure
            % not with control run, but with first run (which includes
            % already the reduced hazard)
            fprintf('SPECIAL CASE FOR BARISAL: Overwrite locations which are not within regional scope of measure with damage from EDS 1!\n')
            EDS(measure_i) = climada_measure_regional_scope(EDS(measure_i),within_scope,EDS(1));
            %EDS(measure_i) = climada_measure_regional_scope(EDS(measure_i),within_scope,EDS(n_measures+1));
        end
    end
        
    if measures.hazard_high_frequency_cutoff(measure_i)>0
        % apply frequency cutoff
        [~,exceedence_freq,~,~,event_index_out] =...
            climada_damage_exceedence(EDS(measure_i).damage,EDS(measure_i).frequency,EDS(measure_i).event_ID);
        cutoff_pos      = exceedence_freq>measures.hazard_high_frequency_cutoff(measure_i);
        event_index_out = event_index_out(cutoff_pos);
        %%fprintf('cutoff %i events\n',length(event_index_out));
        [~,pos] = ismember(event_index_out,EDS(measure_i).event_ID);
        EDS(measure_i).damage(pos) = 0;
        % previous two lines do the same as (but much faster)
        % for index_i=1:length(event_index_out)
        %     EDS_index_i=find(EDS(measure_i).event_ID==event_index_out(index_i));
        %     EDS(measure_i).damage(EDS_index_i)=0;
        % end % index_i
    end
        
    % apply risk transfer
    if measures.risk_transfer_attachement(measure_i)+measures.risk_transfer_cover(measure_i)>0
        damage_in_layer = min(max(EDS(measure_i).damage-measures.risk_transfer_attachement(measure_i),0),...
            measures.risk_transfer_cover(measure_i));
        ED_risk_transfer(measure_i) = full(sum(damage_in_layer.*EDS(measure_i).frequency)); % ED in layer
        EDS(measure_i).damage       = max(EDS(measure_i).damage-damage_in_layer,0); % remaining damage after risk transfer
    end % apply risk transfer

 %end % all except last run
    
    if hazard_switched
        % always switch back, to avoid troubles if hazard passed as struct
        hazard = orig_hazard; % restore
        cprintf([0 0 1],'NOTE: switched hazard back\n');
        if measure_i < n_measures
            orig_hazard=[]; % free up memory
        end
    end
    if assets_switched
        % always switch back, to avoid troubles if hazard passed as struct
        entity.assets = orig_assets; % restore
        fprintf('NOTE: switched assets back\n');
        %orig_assets = []; % free up memory
        assets_switched = 0;
    end
    if entity_switched
        % always switch back, to avoid troubles if entity passed as struct
        entity = orig_entity; % restore
        cprintf([0 0 1],'NOTE: switched entity back\n');
        if measure_i < n_measures
            orig_entity=[]; % free up memory
        end
    end
end % measure_i

measures_impact.EDS = EDS; % store for output

% calculate the cost/benefit ratio also here (so we have all results in one
% -------------------------------------------------------------------------

n_measures = length(measures.cost);
ED         = zeros(1,n_measures); % init
n_years    = climada_global.future_reference_year - climada_global.present_reference_year + 1;

if ~isempty(measures_impact_reference)
    % calculate reference ED (expected damage) no measures
    reference_ED_no_measures = full(sum(measures_impact_reference.EDS(end).damage.*...
        measures_impact_reference.EDS(end).frequency));
end

% get the discound rates for years:
present_year_pos = find(entity.discount.year==climada_global.present_reference_year); % present year
future_year_pos  = find(entity.discount.year==climada_global.future_reference_year); % future year
discount_rates   = entity.discount.discount_rate(present_year_pos:future_year_pos);
ED_no_measures   = full(sum(EDS(end).damage .* EDS(end).frequency)); % calculate annual expected damage
ED(end+1)        = ED_no_measures; % store the ED with no measures at the end (convention)

% time evolution of benefits etc. (linear if climada_global.impact_time_dependence=1)
time_dependence=(0:n_years-1).^climada_global.impact_time_dependence/(n_years-1)^climada_global.impact_time_dependence;

for measure_i = 1:n_measures
    
    % first, calculate the ED (exepected damage) perspective
    ED(measure_i)          = full(sum(EDS(measure_i).damage .* EDS(measure_i).frequency)); % calculate annual expected damage
    ED_benefit(measure_i)  = ED(end) - ED(measure_i);
    
    % store damage frequency curve (DFC), for information only
    DFC(measure_i,:)       = climada_EDS_DFC_report(EDS(measure_i),0,'lean');
    
    % costs are costs as in measures table plus expected damage (for risk transfer only)
    ED_cb_ratio(measure_i) = (measures.cost(measure_i) + ED_risk_transfer(measure_i)) /ED_benefit(measure_i);
    
    % second, calculate the NPV (net present value perspective)
    if isempty(measures_impact_reference)
        % no reference, hence we assume a steady-state, means we see the ED for each year from present to future
        benefits       = ones(1,n_years)*(ED(end)-ED(measure_i)); % same benefit each year
        risk_transfers = ones(1,n_years)*ED_risk_transfer(measure_i); % same risk transfer costs each year
    else
        % time evolution of benefits
        present_benefit = measures_impact_reference.ED(end) - measures_impact_reference.ED(measure_i);
        future_benefit  = ED(end)-ED(measure_i);
        benefits        = present_benefit+(future_benefit-present_benefit)*time_dependence;
        % and similarly for risk transfer costs
        present_risk_transfer = measures_impact_reference.ED_risk_transfer(measure_i);
        future_risk_transfer  = ED_risk_transfer(measure_i);
        risk_transfers        = present_risk_transfer+(future_risk_transfer-present_risk_transfer)*time_dependence;
        
        % old, linear only, kept for backward compatibility, same as climada_global.impact_time_dependence=1
        %         d_benefit       = (future_benefit-present_benefit)/(n_years-1);
        %         d_benefits      = [0,ones(1,n_years-1)*d_benefit];
        %         benefits        = present_benefit+cumsum(d_benefits); % linear increase
        %         % and similarly for risk transfer costs
        %         present_risk_transfer = measures_impact_reference.ED_risk_transfer(measure_i);
        %         future_risk_transfer  = ED_risk_transfer(measure_i);
        %         d_risk_transfer       = (future_risk_transfer-present_risk_transfer)/(n_years-1);
        %         d_risk_transfers      = [0,ones(1,n_years-1)*d_risk_transfer];
        %         risk_transfers        = present_risk_transfer+cumsum(d_risk_transfers);
        
    end
    
    % discount the benefits
    benefit(measure_i)       = climada_NPV(benefits,discount_rates);
    risk_transfer(measure_i) = climada_NPV(risk_transfers,discount_rates); % discount the risk transfer costs
    % costs are costs as in measures table plus expected damage (for risk transfer only)
    cb_ratio(measure_i)      = (measures.cost(measure_i)+risk_transfer(measure_i))/benefit(measure_i);
    
end % measure_i

% calculate the NPV of the full unaverted damages, too
% TCR stands for total climate risk
if isempty(measures_impact_reference)
    NPV_total_climate_risk = climada_NPV(ones(1,n_years)*ED(end), discount_rates);
else
    % time evolution of risk
    present_TCR = measures_impact_reference.ED(end);
    future_TCR  = ED(end);
    TCRs        = present_TCR+(future_TCR-present_TCR)*time_dependence;
    % old, linear only, kept for backward compatibility, same as climada_global.impact_time_dependence=1
    % d_TCR       = (future_TCR-present_TCR)/(n_years-1);
    % d_TCRs      = [0,ones(1,n_years-1)*d_TCR];
    % TCRs        = present_TCR+cumsum(d_TCRs); % linear increase
    NPV_total_climate_risk = climada_NPV(TCRs, discount_rates);
end

% store in measures
measures_impact.ED               = ED;
measures_impact.DFC              = DFC; % info only
measures_impact.ED_benefit       = benefit;
measures_impact.ED_risk_transfer = ED_risk_transfer;
measures_impact.ED_cb_ratio      = cb_ratio;
measures_impact.benefit          = benefit;
measures_impact.risk_transfer    = risk_transfer;
measures_impact.cb_ratio         = cb_ratio;
measures_impact.NPV_total_climate_risk = NPV_total_climate_risk;
measures_impact.peril_ID         = hazard.peril_ID;
measures_impact.measures         = measures; % store measures into res, so we have a complete set

% prepare annotation
[~,hazard_name]   = fileparts(EDS(1).hazard.filename);
[~,assets_name]   = fileparts(EDS(1).assets.filename);

if ~isfield(measures,'filename')
    measures_name = ['Measures_' datestr(now,'yymmddHHMM')];
else
    [~,measures_name] = fileparts(measures.filename);
end
if strcmp(measures_name,assets_name),measures_name='m';end
measures_impact.title_str = sprintf('%s @ %s | %s',measures_name,assets_name,hazard_name);

save_filename = strrep(measures_impact.title_str,' ',''); % for filename
save_filename = strrep(save_filename,'_',''); % for filename
save_filename = strrep(save_filename,'@','_'); % for filename
save_filename = strrep(save_filename,'|','_'); % for filename

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

save_filename = [module_data_dir filesep 'results' filesep save_filename];

measures_impact.filename = save_filename;

% last but not least, calculate risk premium
measures_impact.risk_premium_fgu=measures_impact.NPV_total_climate_risk/measures_impact.EDS(1).Value;
fprintf('total climate risk premium (fgu - only a proxy):        %f%%\n',...
    measures_impact.risk_premium_fgu*100);
measures_impact.risk_premium_net=(measures_impact.NPV_total_climate_risk-...
    sum(measures_impact.ED_benefit(measures_impact.cb_ratio<1)))/measures_impact.EDS(1).Value;
fprintf('total climate risk premium (net of effective measures): %f%%\n',...
    measures_impact.risk_premium_net*100);
fprintf('total climate risk premium reduction (relative)        %2.2f%%\n',...
    -(measures_impact.risk_premium_net/measures_impact.risk_premium_fgu-1)*100);
measures_impact.risk_premium_comment='total climate risk divided by the total value of assets';

% save(save_filename,'measures_impact')
% fprintf('results written to %s\n',save_filename);

if map_risk_premium
    % plot the total climate risk premium
    
    % risk premium fgu at each centroid
    centroid_risk_premium_fgu=measures_impact.EDS(end).ED_at_centroid;
    
    % risk premium net of effective measures at each centroid
    d_risk_premium_net=centroid_risk_premium_fgu*0; % init
    effective_measures_pos=find(measures_impact.cb_ratio<1);
    for i=1:length(effective_measures_pos) % sum over all benefits with c/b <1
        d_risk_premium_net=d_risk_premium_net+(centroid_risk_premium_fgu-...
            measures_impact.EDS(effective_measures_pos(i)).ED_at_centroid);
    end % i
    centroid_risk_premium_net=centroid_risk_premium_fgu-d_risk_premium_net;
    % divide by total Value ast each centroid (to obtain risk premium)
    centroid_risk_premium_fgu=centroid_risk_premium_fgu./measures_impact.EDS(end).assets.Value*100;
    centroid_risk_premium_net=centroid_risk_premium_net./measures_impact.EDS(end).assets.Value*100;
    
    % define a common color axis
    caxis_range=[min(min(centroid_risk_premium_fgu),min(centroid_risk_premium_net)) ...
        max(max(centroid_risk_premium_fgu),max(centroid_risk_premium_net))];
    
    % plot on two panes, figure whether wider extent in lon or lat
    dlon=max(measures_impact.EDS(end).assets.lon)-min(measures_impact.EDS(end).assets.lon);
    dlat=max(measures_impact.EDS(end).assets.lat)-min(measures_impact.EDS(end).assets.lat);
    if dlat>dlon,subplot(1,3,1),else subplot(3,1,1);end
    climada_color_plot(centroid_risk_premium_fgu,...
        measures_impact.EDS(end).assets.lon,...
        measures_impact.EDS(end).assets.lat,'none','risk premium fgu [%]',...
        'pcolor','linear',199,1,caxis_range);
    if dlat>dlon,subplot(1,3,2),else subplot(3,1,2);end
    climada_color_plot(centroid_risk_premium_net,...
        measures_impact.EDS(end).assets.lon,...
        measures_impact.EDS(end).assets.lat,'none','risk premium net [%]',...
        'pcolor','linear',199,1,caxis_range);
    set(gcf,'Color',[1 1 1]);
    if dlat>dlon,subplot(1,3,3),else subplot(3,1,3);end
    climada_color_plot(-(centroid_risk_premium_net./centroid_risk_premium_fgu-1)*100,...
        measures_impact.EDS(end).assets.lon,...
        measures_impact.EDS(end).assets.lat,'none','risk premium reduction [rel %]',...
        'pcolor','linear',199,1);
    set(gcf,'Color',[1 1 1]);
    
end % map_risk_premium

return
