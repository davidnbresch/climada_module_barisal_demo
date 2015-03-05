

%% create barisal BCC entity


global climada_global
if ~climada_init_vars,return;end % init/import global variables

% read BCC boundaries
GIS_dir  = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\Barisal_GIS\WGS1984';
BCC_file = 'CityCorporationAreaPolyBCC_P.shp'; 
BCC      = climada_shaperead([GIS_dir filesep BCC_file],0,1,0,1); % 0 for no-save

% load barisal division entity
barisal_entity_mat = 'barisal_01x01.mat';
load([climada_global.data_dir filesep 'entities' filesep barisal_entity_mat]);
figure
climada_plot_entity_assets(entity)
% climada_cut_out_GDP_entity
hold on 
% figure
h(1) = plot(BCC(1).X,BCC(1).Y,'color',[240 128 128]/255);

% find asset points within the boundary
polygon = [BCC.X' BCC.Y'];
polygon(isnan(polygon(:,1)),:) = [];
cn = inpoly([entity.assets.lon' entity.assets.lat'],polygon);
fprintf('--> %d from %d asset points are within the BCC boundaries. \n', sum(cn), numel(cn))

% cut out specific points
entity_ori = entity;
entity.assets.lon           = entity.assets.lon(cn);
entity.assets.lat           = entity.assets.lat(cn);
entity.assets.Value         = entity.assets.Value(cn);
entity.assets.Deductible    = entity.assets.Deductible(cn);
entity.assets.Cover         = entity.assets.Cover(cn);
entity.assets.DamageFunID   = entity.assets.DamageFunID(cn);
        


figure
climada_plot_entity_assets(entity)
hold on
plot(BCC(1).X,BCC(1).Y,'color',[240 128 128]/255);

% save barisal BCC entity
save([climada_global.modules_dir filesep 'barisal_demo' filesep 'data' filesep 'entities' filesep 'Barisal_BCC_1km.mat'], 'entity')

% scale to 100 total asset value
tot_val = sum(entity.assets.Value);
entity.assets.Value      = entity.assets.Value/tot_val*100;
entity.assets.Deductible = entity.assets.Deductible/tot_val*100;
entity.assets.Cover      = entity.assets.Cover/tot_val*100;
save([climada_global.modules_dir filesep 'barisal_demo' filesep 'data' filesep 'entities' filesep 'Barisal_BCC_1km_100.mat'], 'entity')


%% create centroids for wind

% 80 km grid for bangaldesh
resolution_km  = 80;
centroids_rect = [87 93 19 26];
centroids_grid = climada_generate_HR_centroids(centroids_rect, resolution_km);

% 20 km grid for barisal district
resolution_km2   = 20;
centroids_rect2 = [89.88-0.18*2 91.2 21.88-0.18*2 23.5];
centroids_grid2 = climada_generate_HR_centroids(centroids_rect2, resolution_km2);

% merge centroids
[a b] = max(size(entity.assets.lon));
centroids.lon = cat(b, entity.assets.lon, centroids_grid2.lon, centroids_grid.lon);
centroids.lat = cat(b, entity.assets.lat, centroids_grid2.lat, centroids_grid.lat);
centroids.centroid_ID = 1:numel(centroids.lon);
centroids.comment = ['handmade, ' datestr(now)];

% save centroids
save([climada_global.modules_dir filesep 'barisal_demo' filesep 'data' filesep 'system' filesep 'Barisal_BCC_centroids.mat'], 'centroids')

% to visually check the focus region
close all
fig = climada_figuresize(0.5,0.7);
climada_plot_world_borders
axis(centroids_rect)
axis equal
axis(centroids_rect)

hold on
plot(centroids_grid2.lon, centroids_grid2.lat,'+g')
plot(centroids_grid.lon, centroids_grid.lat,'+r')
plot(BCC(1).X,BCC(1).Y,'color',[240 128 128]/255);
plot(BCC(1).X(1),BCC(1).Y(1),'+','color',[240 128 128]/255,'markersize',10, 'linewidth', 3);
title(sprintf('Barisal, %d centroids', numel(centroids.lon)))


% print figure as pdf
foldername = [climada_global.modules_dir filesep 'barisal_demo' filesep 'data' filesep 'system' filesep 'Barisal_centroids.pdf'];
print(fig,'-dpdf',foldername)



% alternative
% [centroids_new,entity_new,polygon] = climada_cut_out_GDP_entity(entity,centroids,polygon);



