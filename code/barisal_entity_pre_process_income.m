function entity = barisal_entity_pre_process_income(entity)

%% create income column in entity.assets structure, so that the income 
% corresponds to the residential buildings
% - 4 different income figures for Juphri, Katcha, Pucca and Semi-Pucca houses
% - Pucca and Semi-Pucca house values are distributed to normal and +30cm houses, 
%   whereas we need to assign the income value relative to the portion of normal or +30cm houses



% unique asset categories
categories_uni = unique(entity.assets.Category);

% residential categories that we have income information for
categories_residential = {'buildings_Juphri' 'buildings_Katcha' 'buildings_Pucca' 'buildings_Semi_Pucca'};

% income field names
flds_names_income = {'Income_of_juphri_residents' 'Income_of_katcha_residents' 'Income_of_semipucca_residents' 'Income_of_pucca_residents'};

% all categories that contain residential categories
indx_residential = strfind(categories_uni,'Residential');
indx_residential = find(~cellfun(@isempty,indx_residential));

% all categories that are part of the +30 cm category
indx_30 = strfind(categories_uni,'30');
indx_30 = find(~cellfun(@isempty,indx_30));

% find asset categories where we have income information for
for c_i = 1:numel(categories_residential)
    indx  = strfind(categories_uni,categories_residential{c_i});
    indx  = find(~cellfun(@isempty,indx));
    indx  = indx(ismember(indx, indx_residential));  
    %indx  = indx(~ismember(indx, indx_30));  
    %fprintf('%s, %s\n',categories_residential{c_i},categories_uni{indx});
    income2category_indx{c_i} = indx;
end

% init income vector
entity.assets.income = zeros(size(entity.assets.lon));

for c_i = 1:numel(categories_residential)
    
    Value_temp = [];
    % find all values that belong to this income information
    indx = income2category_indx{c_i};
    for ii = 1:numel(indx)
        assets_indx = strcmp(entity.assets.Category,categories_uni{indx(ii)});
        Value_temp(:,ii) = entity.assets.Value(assets_indx);
    end
    % add up all values for this income information (e.g. residential pucca & residential pucca +30cm)
    Value_temp_tot = sum(Value_temp,2);
    
    % assign income value relative to building values (portion)
    for ii = 1:numel(indx)
        assets_indx = strcmp(entity.assets.Category,categories_uni{indx(ii)});
        portion = Value_temp(:,ii)./Value_temp_tot;
        portion(isnan(portion)) = 0;
        entity.assets.income(assets_indx) = entity.assets.(flds_names_income{c_i})(1:sum(assets_indx)) .*portion;
    end
    
end





