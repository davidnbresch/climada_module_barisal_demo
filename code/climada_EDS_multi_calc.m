function [EDS, hazards, entities] = climada_EDS_multi_calc(entities,hazards,EDS_save_file,sync_check,check_plots)

% damage calc for combination of hazards and entities

global climada_global
if ~climada_init_vars, return; end

if ~exist('entities',       'var'), entities        = '';   end
if ~exist('hazards',        'var'), hazards         = '';   end
if ~exist('EDS_save_file',  'var'), EDS_save_file   = '';   end
if ~exist('sync_check',     'var'), sync_check      = 1;    end
if ~exist('check_plots',    'var'), check_plots 	= 0;    end

module_data_dir=[fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];

% prompt for entity
if isempty(entities)
    entities=[module_data_dir filesep 'entities' filesep '*.mat'];
    [fN, fP] = uigetfile(entities, 'Select entities:','MultiSelect','on');
else
    if ischar(entities) % only one filename passed
        [fP, fN] = fileparts(entities);
        if ~climada_check_matfile(entities)
            climada_entity_read(entities,'NOENCODE')
        end
        fN = [fN '.mat'];
    elseif iscellstr(entities) % multiple filenames passed in 
        for e_i = 1:length(entities)
            [fP, fN{e_i}] = fileparts(entities{e_i});
            fN{e_i} = [fN{e_i} '.mat'];
        end   
    end
end
clear entities

if isequal(fN,0) || isequal(fP,0)
    return; % cancel
else
    if iscell(fN)
        fmt_str = '%s';
        for e_i = 1 : length(fN)
            msg_str = sprintf('loading entity file: %s',fN{e_i});
            fprintf(fmt_str,msg_str)
            fmt_str = [repmat('\b',size(msg_str)) '%s'];
            load(fullfile(fP,fN{e_i}));
            if ~exist('entity','var'), cprintf([1 0 0],'ERROR: invalid entity file \n'); return; end
            entities(e_i) = entity;
            % for consistency
            entities(e_i).assets.filename = fullfile(fP,fN{e_i});
        end
        fprintf(fmt_str,'entity files loaded \n')
    else
        fprintf('loading entity file: %s... ',fN);
        load(fullfile(fP,fN));
        if ~exist('entity','var'), cprintf([1 0 0],'ERROR: invalid entity file \n'); return; end
        entities = entity;
        % for consistency
        entities(e_i).assets.filename = fullfile(fP,fN{e_i});
        fprintf('done \n')
    end
    clear entity
end


for e_i = 1:length(entities)
    if ~isfield(entities(e_i).assets,'reference_year')
        entities(e_i).assets.reference_year = climada_global.present_reference_year;
    end
    
    flds = fieldnames(entities(e_i).assets);
    rm_flds = {};
    for fld_i = 1:length(flds)
        if strfind(flds{fld_i},'Value_')
            rm_flds(end+1) = flds(fld_i);
            
            tmp_assets          = entities(e_i).assets;
            %             tmp_assets          = rmfield(entity(e_i).assets,flds{fld_i});
            
            tmp_assets.Value    = entities(e_i).assets.(flds{fld_i});
            C = strsplit(flds{fld_i},'_');
            tmp_assets.reference_year = str2double(C{2});
            
            entities(end+1).damagefunctions   = entities(e_i).damagefunctions;
            entities(end).measures            = entities(e_i).measures;
            entities(end).discount            = entities(e_i).discount;
            entities(end).assets              = tmp_assets;
            
        end
        clear tmp_assets
    end
end

for e_i = 1:length(entities)
    a_flds = fieldnames(entities(e_i).assets);
    for a_fld_i = 1:length(a_flds)
        if ismember(a_flds{a_fld_i},rm_flds)
            entities(e_i).assets  = rmfield(entities(e_i).assets,a_flds{a_fld_i});
        end
    end
end

clear e_i a_flds a_fld_i flds fld_i rm_flds

% fields necessary for damage calc

% prompt for hazard set
if isempty(hazards)
    hazards=[module_data_dir filesep 'hazards' filesep '*.mat'];
    [fN, fP] = uigetfile(hazards, 'Select hazard:','MultiSelect','on');
else
    if ischar(hazards) % only one filename passed
        [fP, fN] = fileparts(hazards);
        fN = [fN '.mat'];
    elseif iscellstr(hazards) % multiple filenames passed in 
        for h_i = 1:length(hazards)
            [fP, fN{h_i}] = fileparts(hazards{h_i});
            fN{h_i} = [fN{h_i} '.mat'];
        end
    end
end
clear hazards

if isequal(fN,0) || isequal(fP,0)
    return; % cancel
else
    if iscell(fN)
        fmt_str = '%s'; hazards = {};
        for h_i = 1 : length(fN)
            msg_str = sprintf('loading hazard file: %s',fN{h_i});
            fprintf(fmt_str,msg_str)
            fmt_str = [repmat('\b',size(msg_str)) '%s'];
            load(fullfile(fP,fN{h_i}));
            if ~exist('hazard','var'), cprintf([1 0 0],'ERROR: invalid hazard file \n'); return; end
            
            if ~isfield(hazard,'reference_year')
                cprintf([1 0.5 0],sprintf('WARNING: no reference year for %s - using default %4.0f\n',...
                    fN{h_i},climada_global.present_reference_year));
                hazard.reference_year = climada_global.present_reference_year;
            end
            
            % for consistency
            hazard.filename = fullfile(fP,fN{h_i});
            hazard.comment = strrep(strtok(fN{h_i},'.'),'_',' ');
            hazards{h_i}=hazard;
            
        end
        fprintf(fmt_str,'hazard files loaded \n');
    else
        fprintf('loading hazard file: %s... ',fN);
        load(fullfile(fP,fN));
        if ~exist('hazard','var'), cprintf([1 0 0],'ERROR: invalid hazard file \n'); return; end
        
        if ~isfield(hazard,'reference_year')
            cprintf([1 0.5 0],sprintf('WARNING: no reference year for %s - using default %4.0f\n',...
                fN{h_i},climada_global.present_reference_year));
            hazard.reference_year = climada_global.present_reference_year;
        end
        
        % for consistency
        hazard.filename = fullfile(fP,fN);
        hazard.comment = strrep(strtok(fN{h_i},'.'),'_',' ');
        hazards={hazard};
        fprintf('done \n')
    end
    clear hazard
end

ed_i = 0;

for e_i = 1: length(entities)
    for h_i = 1: length(hazards)
        hazard_i = hazards{h_i};
        
        if sync_check && (hazard_i.reference_year ~= entities(e_i).assets.reference_year)
            continue;
        end
        
        ed_i = ed_i +1;
%         fprintf('hazard.comment = ''%s'' \n- would you like to use this for the EDS annotation?\n',hazard(h_i).comment)
%         h_str = input('to continue press ''y'', otherwise type a new hazard title:\n','s');
%
%         if strcmp(h_str,'y'),   h_str = hazard(h_i).comment;        end
%
%        	fprintf('entity.assets.comment = ''%s'' \n- would you like to use this for the EDS annotation?\n',entity(e_i).assets.comment)
%         e_str = input('to continue press ''y'', otherwise type a new entity title:\n','s');
%
%         if strcmp(e_str,'y'),   e_str = entity(e_i).assets.comment; end
        
        h_str = hazard_i.comment;
        e_str = entities(e_i).assets.comment;
        
%         annotation  = sprintf('%s | %s (%d)',h_str,e_str,entities(e_i).assets.reference_year);
        annotation  = sprintf(h_str);
        EDS(ed_i)   = climada_EDS_calc(entities(e_i),hazard_i,annotation);
        if EDS(ed_i).ED == 0
            msg = sprintf('WARNING: expected damage equals zero for entity %s and hazard %s, removing from EDS structure\n',...
                e_str,h_str);
            cprintf([1 0.5 0],msg)
            EDS(ed_i) = [];
            ed_i = ed_i - 1;
        end
    end
end

if ~isempty(EDS_save_file) && ischar(EDS_save_file) && ~strcmp(EDS_save_file,'NO_SAVE')
    save(EDS_save_file,'EDS')
elseif isempty(EDS_save_file)
    uisave('EDS',[module_data_dir filesep 'results' filesep 'multi_EDS.mat']);
end

if check_plots
    figure('color','w','name',annotation)
    climada_ED_plot(EDS,0);
end

