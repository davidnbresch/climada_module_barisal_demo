function entity = barisal_entity_pre_process(entity_in)
% barisal_entity_pre_process
% MODULE:
%   barisal_demo
% NAME:
%   barisal_entity_pre_process
% PURPOSE:
%   Pre process entity into climada entity format from input as provided
%   by Barisal/Ecorys team, where assetinformation is in a matrix 
%   (lon, lat, Ward Nr, and two columns for every category (one with 
%   category name ASSETS and second with category name DAMAGE FUNCTION).
%   check that column headers have unique names, otherwise data is lost.
% CALLING SEQUENCE:
%   entity = barisal_entity_pre_process(entity_in)
% EXAMPLE:
%   entity = barisal_entity_pre_process(entity_in)
% INPUTS:
%   entity_in: entity as read from Ecorys input in matrix format
%       climada_xlsread('no',ecorys_entity_file_xls,'Floods_2014',1);
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   entity: a climada entity structure
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150909, add documentation
%-

% entity_in = entity;
entity = entity_in;

entity.assets.Category = [];

% find fields with asset information (for all the different asset categories)
flds = fieldnames(entity_in.assets);
ndx_asset_flds  = strfind(flds,'ASSETS');
ndx_asset_flds  = find(~cellfun(@isempty,ndx_asset_flds));

for asset_fld_i = 1:length(ndx_asset_flds)
    
    % copy asset VALUE to temporary field
    value_tmp  = entity_in.assets.(flds{ndx_asset_flds(asset_fld_i)});
    if iscell(value_tmp)   
        ndx = cellfun(@ischar,value_tmp);
        value_tmp(ndx) ={0};
        value_tmp = cell2mat(value_tmp);
    end
    % copy existing VALUE
    value_part = entity.assets.Value;
    % add up existing VALUE and temp VALUE
    if all(isnan(value_part)) 
        % special case for first time
        value = value_tmp;
    else 
        value = [value_part; value_tmp];
    end
    % rewrite total to VALUE field
    entity.assets.Value = value;
    
    
    % copy DamageFunID to temporary field
    DamageFunID_tmp = entity_in.assets.(flds{ndx_asset_flds(asset_fld_i)+1});
    if iscell(DamageFunID_tmp)   
        ndx = cellfun(@ischar,DamageFunID_tmp);
        DamageFunID_tmp(ndx) ={0};
        DamageFunID_tmp = cell2mat(DamageFunID_tmp);
    end
    
    % copy existing DamageFunID
    DamageFunID_part = entity.assets.DamageFunID;
    % add up existing VALUE and temp VALUE
    if all(isnan(DamageFunID_part)) 
        % special case for first time
        DamageFunID = DamageFunID_tmp;
    else 
        DamageFunID = [DamageFunID_part; DamageFunID_tmp];
    end
    % rewrite total to DamageFunID field
    entity.assets.DamageFunID = DamageFunID;
    
    % add Category field
    category_names = {strrep(flds{ndx_asset_flds(asset_fld_i)},', ASSETS','')};
    if isfield(entity_in.assets,'X')
        category_names = repmat(category_names,length(entity_in.assets.X),1);
    else
        category_names = repmat(category_names,length(entity_in.assets.lon),1);
    end
    
    % copy existing categories
    category_part = entity.assets.Category;
    if isempty(category_part)
        category      = category_names;
    else
        category      = {category_part{:} category_names{:}};
    end
    entity.assets.Category = category;

    fprintf('\t - added %s\n', flds{ndx_asset_flds(asset_fld_i)})
end


% find fields to be multiplied
invalid_fld_1 = strfind(flds,'ASSETS');
invalid_fld_2 = strfind(flds,'DAMAGE FUNCTION');
invalid_fld_1 = ~cellfun(@isempty,invalid_fld_1);
invalid_fld_2 = ~cellfun(@isempty,invalid_fld_2);
invalid_flds  = logical(invalid_fld_1+invalid_fld_2);

% special cases for fieldnames that do not need to be multiplied
ndx_flds_to_copy = ~(invalid_flds+strcmp(flds,'filename')+strcmp(flds,'hazard')+strcmp(flds,'Value')+strcmp(flds,'DamageFunID')+strcmp(flds,'Category')+strcmp(flds,'Income'));
ndx_flds_to_copy = find(ndx_flds_to_copy);

for fld_i = 1:length(ndx_flds_to_copy)
    val_temp  = entity_in.assets.(flds{ndx_flds_to_copy(fld_i)});
    val_total = repmat(val_temp,length(ndx_asset_flds),1);
    entity.assets.(flds{ndx_flds_to_copy(fld_i)}) = val_total;
end
    

% finally remove both fields (assets and DMG function)
for asset_fld_i = length(ndx_asset_flds):-1:1
    if isfield(entity.assets,flds{ndx_asset_flds(asset_fld_i)+1})
        entity.assets = rmfield(entity.assets,flds{ndx_asset_flds(asset_fld_i)+1});
    end
    entity.assets = rmfield(entity.assets,flds{ndx_asset_flds(asset_fld_i)});
end  
    
    
    
    