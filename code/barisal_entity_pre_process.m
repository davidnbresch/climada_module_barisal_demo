function entity_out = barisal_entity_pre_process(entity_in)

% entity_in = entity;
entity_out = entity_in;

entity_out.assets.Category = [];

% find fields with asset information (for all the different asset categories)
names = fieldnames(entity_in.assets);
indx_asset_fields  = strfind(names,'ASSETS');
indx_asset_fields  = find(~cellfun(@isempty,indx_asset_fields));


for value_i = 1:length(indx_asset_fields)
    
    % copy asset VALUE to temporary field
    value_temp  = getfield(entity_in.assets, names{indx_asset_fields(value_i)});
    if iscell(value_temp)   
        indx = cellfun(@ischar,value_temp);
        value_temp(indx) ={0};
        value_temp = cell2mat(value_temp);
    end
    % copy existing VALUE
    value_part = getfield(entity_out.assets, 'Value');
    % add up existing VALUE and temp VALUE
    if all(isnan(value_part)) 
        % special case for first time
        value = value_temp;
    else 
        value = [value_part; value_temp];
    end
    % rewrite total to VALUE field
    entity_out.assets = setfield(entity_out.assets,'Value',value);
    
    
    % copy DamageFunID to temporary field
    DamageFunID_temp = getfield(entity_in.assets, names{indx_asset_fields(value_i)+1});
    if iscell(DamageFunID_temp)   
        indx = cellfun(@ischar,DamageFunID_temp);
        DamageFunID_temp(indx) ={0};
        DamageFunID_temp = cell2mat(DamageFunID_temp);
    end
    
    % copy existing DamageFunID
    DamageFunID_part = getfield(entity_out.assets, 'DamageFunID');
    % add up existing VALUE and temp VALUE
    if all(isnan(DamageFunID_part)) 
        % special case for first time
        DamageFunID = DamageFunID_temp;
    else 
        DamageFunID = [DamageFunID_part; DamageFunID_temp];
    end
    % rewrite total to DamageFunID field
    entity_out.assets = setfield(entity_out.assets,'DamageFunID',DamageFunID);
    
    % add Category field
    category_names = {strrep(names{indx_asset_fields(value_i)},'_ASSETS','')};
    if isfield(entity_in.assets,'X')
        category_names = repmat(category_names,length(entity_in.assets.X),1);
    else
        category_names = repmat(category_names,length(entity_in.assets.lon),1);
    end
    
    % copy existing categories
    category_part = getfield(entity_out.assets, 'Category');
    if isempty(category_part)
        category      = category_names;
    else
        category      = {category_part{:} category_names{:}};
    end
    entity_out.assets = setfield(entity_out.assets,'Category',category');

    fprintf('\t - added %s\n', names{indx_asset_fields(value_i)})
end


% find fields to be multiplied
invalid_field_1 = strfind(names,'ASSETS');
invalid_field_2 = strfind(names,'DMG_Function');
invalid_field_1 = ~cellfun(@isempty,invalid_field_1);
invalid_field_2 = ~cellfun(@isempty,invalid_field_2);
invalid_fields = logical(invalid_field_1+invalid_field_2);

% special cases for fieldnames that do not need to be multiplied
indx_to_be_multiplied_fields = ~(invalid_fields+strcmp(names,'filename')+strcmp(names,'hazard')+strcmp(names,'Value')+strcmp(names,'DamageFunID')+strcmp(names,'Category'));
indx_to_be_multiplied_fields = find(indx_to_be_multiplied_fields);

for field_i = 1:length(indx_to_be_multiplied_fields)
    val_temp  = getfield(entity_in.assets, names{indx_to_be_multiplied_fields(field_i)});
    val_total = repmat(val_temp,length(indx_asset_fields),1);
    entity_out.assets = setfield(entity_out.assets,names{indx_to_be_multiplied_fields(field_i)},val_total);
end
    

% finally remove both fields (assets and DMG function)
for value_i = length(indx_asset_fields):-1:1
    entity_out.assets = rmfield(entity_out.assets,names{indx_asset_fields(value_i)+1});
    entity_out.assets = rmfield(entity_out.assets,names{indx_asset_fields(value_i)});
end  
    
    
    
    