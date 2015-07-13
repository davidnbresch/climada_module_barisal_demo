function entity = barisal_ED_find_most_vulnerable(entity, EDS, criterion_A, criterion_B)

%% compare AED with income

% set to small income numbers to zero
entity.assets.income(entity.assets.income<1) = 0;

% calculate relative AED, relative to income
AED_rel = EDS(1).ED_at_centroid(:) ./entity.assets.income(:) ;
% AED_rel = measures_impact3(1).EDS(1).ED_at_centroid(:) ./entity.assets.income(:) ;
AED_rel(isinf(AED_rel)) = 0;
AED_rel(isnan(AED_rel)) = 0;

% find top-most relative AEDs
% Y = prctile(AED_rel,75) 
[AED_sort AED_sort_indx] = sort(AED_rel,'descend');

% number of centroids with residential categories
n_residential = 6*7753;

% find most vulnerable buildings
criterion_absolute_A      = int64(criterion_A*n_residential);
criterion_absolute_B      = int64(criterion_B*n_residential);
indx_most_vuln_building_A = AED_sort_indx(1:criterion_absolute_A);
indx_most_vuln_building_B = AED_sort_indx(1:criterion_absolute_B);

% index with centroids that point to the most vulnerable residential buildings
indx_all = 1:numel(entity.assets.lon);
indx_most_vuln_building_A = ismember(indx_all, indx_most_vuln_building_A);
indx_most_vuln_building_B = ismember(indx_all, indx_most_vuln_building_B);

% sum(indx_most_vuln_building)

% load resilient_buildings_zone_B
barisal_data_dir= [fileparts(fileparts(mfilename('fullpath'))) filesep 'data'];
load([barisal_data_dir filesep 'entities' filesep 'Measures_package' filesep 'resilient_buildings_zones_B.mat'])

% find lat/lon in zone B
indx_zone_B = inpoly([entity.assets.lon entity.assets.lat],[resilient_buildings_zone_B.lon' resilient_buildings_zone_B.lat']);

% create vector with 
%     - zeros (not most vulnerable)
%     - A (most vulnerable buildings in zone A) and 
%     - B (most vulnerable buildings in zone B)
most_vuln_zone_B   = logical(indx_most_vuln_building_B .* indx_zone_B');
most_vuln_zone_A   = logical(indx_most_vuln_building_A .* ~indx_zone_B');
most_vuln_building = cell(size(entity.assets.lon));
most_vuln_building(most_vuln_zone_A) = repmat({'A'},sum(most_vuln_zone_A),1);
most_vuln_building(most_vuln_zone_B) = repmat({'B'},sum(most_vuln_zone_B),1);
entity.assets.most_vuln_building     = most_vuln_building;

climada_figuresize(0.4,0.4)
plot(entity.assets.lon(most_vuln_zone_A), entity.assets.lat(most_vuln_zone_A),'or','markersize',2)
hold on
plot(entity.assets.lon(most_vuln_zone_B), entity.assets.lat(most_vuln_zone_B),'ob','markersize',2)
title(sprintf('%d%% most vulnerable res. buildings in Zone A (red) and %d%% in Zone B (blue)',criterion_A*100,criterion_B*100),'fontsize',10)
axis equal

%% extend to other buildings

% % unique asset categories
% [categories_uni ia ic] = unique(entity.assets.Category);

% residential categories that we have income information for
categories_residential = {'Residential_buildings_Pucca_ASSETS'...
                          'Residential_buildings_Semi_Pucca_ASSETS'...
                          'Residential_buildings_Katcha_ASSETS'...
                          'Residential_buildings_Pucca_ASSETS_30_cm_elevation_'...
                          'Residential_buildings_Semi_Pucca_ASSETS_30_cm_elevation_'};

more_building_categories = {'Commercial' 'industry' 'public'};

% loop over construction types
for c_i = 1:numel(categories_residential)
    
    indx = strcmp(entity.assets.Category, categories_residential{c_i});
    %sum(indx)
    
    % loop over building categories
    for c_ii = 1:numel(more_building_categories)
        category_name = strrep(categories_residential{c_i},'Residential',more_building_categories{c_ii});
        indx_2        = strcmp(entity.assets.Category, category_name);
        %sum(indx_2)
        entity.assets.most_vuln_building(indx_2) = entity.assets.most_vuln_building(indx);
    end
end
                          



% climada_figuresize(0.4,0.4)
% indxA = strcmp(entity.assets.most_vuln_building,'B') & strcmp(entity.assets.Category, 'public_buildings_Pucca_ASSETS')';
% plot(entity.assets.lon(indxA), entity.assets.lat(indxA),'or','markersize',2)
% hold on
% plot(entity.assets.lon(strcmp(entity.assets.most_vuln_building,'B')), entity.assets.lat(strcmp(entity.assets.most_vuln_building,'B')),'ob','markersize',2)
% title(sprintf('%d%% most vulnerable res. buildings in Zone A (red) and B (blue)',criterion*100),'fontsize',10)
% axis equal






