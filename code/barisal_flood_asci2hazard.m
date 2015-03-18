function hazard = barisal_flood_asci2hazard
% barisal_flood_asci2hazard
% MODULE:
%   barisal_demo
% NAME:
%   barisal_flood_asci2hazard
% PURPOSE:
%   read flood asci files and transform to climada hazard structure, asci
%   files need to be saved in climada/data/hazard/Flood_Barisal/today/...,
%   years given 1983 to 2011 (29 years, 29 events), save hazard in
%   data/hazard
% CALLING SEQUENCE:
%   hazard=barisal_flood_asci2hazard
% EXAMPLE:
%   hazard=barisal_flood_asci2hazard
% INPUTS: 
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS: hazard
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150313, init
%-

hazard = []; %init

global climada_global
if ~climada_init_vars,return;end % init/import global variables


%% read flood asci file
% climada_global.data_dir = '\\CHRB1065.CORP.GWPNET.COM\homes\X\S3BXXW\Documents\lea\climada_git\climada_data';
ascifile = [climada_global.data_dir filesep 'hazards' filesep 'Flood_Barisal' filesep 'today' filesep 'MaxInundationDepths1983.asc'];
delimiter = '\t';
R1 = 6;
C1 = 0;
flood_grid = flipud(dlmread(ascifile,delimiter,R1,C1));
nodata_value = -9999;
flood_grid(flood_grid==nodata_value) = 0;

%calculate lat lon corners
xllcorner = 530285.438;
yllcorner = 504276.750;
cellsize  = 100.0;

% transform for lat lon including shift
[no_row, no_col]   = size(flood_grid);
[lat_min, lon_min] = utm2ll_shift(xllcorner, yllcorner);
[lat_max, lon_max] = utm2ll_shift(xllcorner+cellsize*no_col, yllcorner+cellsize*no_row);

% original conversion from UTM to lat lon
% [lat_min, lon_min] = btm2ll(xllcorner, yllcorner);
% [lat_max, lon_max] = btm2ll(xllcorner+cellsize*no_col, yllcorner+cellsize*no_row);

% create meshgrid
[X, Y ] = meshgrid(linspace(lon_min,lon_max,no_col), ...
                    linspace(lat_min,lat_max,no_row));
                
%% test figure                
% figure
% contourf(X,Y,flood_grid)
% hold on
% climada_plot_world_borders('', '', '', 1);
% cmap = climada_colormap('FL');
% colormap(cmap)
% colorbar
% % caxis([])


%% transform to hazard-structure
% load TC example and use as template for flood hazard
load([climada_global.modules_dir filesep 'barisal_demo' filesep 'data' filesep 'hazards' filesep 'Barisal_BCC_hazard_TC.mat'])
hazard_TC = hazard;

% overwrite template hazard with flood information
no_event   = 29;
hazard.lon = reshape(X,1,no_col*no_row);
hazard.lat = reshape(Y,1,no_col*no_row);  
hazard.centroid_ID = 1:numel(hazard.lon);
hazard.orig_years       = no_event;
hazard.orig_event_count = no_event;
hazard.event_count      = no_event;
hazard.event_ID         = 1:no_event;
hazard.orig_event_flag  = ones(1,no_event);
hazard.yyyy = ones(1,no_event);
hazard.mm   = ones(1,no_event);
hazard.dd   = ones(1,no_event);
hazard.intensity = sparse(no_event,numel(hazard.lon));
hazard.name      = cell(1,no_event);
hazard.frequency = ones(1,no_event)/no_event;
hazard.peril_ID  = 'FL';
hazard.comment   = 'modelled by W+B';
hazard.date      = datestr(now);
hazard.units     = 'm';
hazard.orig_yearset = [];

filename   = 'MaxInundationDepths';
year_start = 1983;

% transform flood heights to hazard intensity
for e_i = 1:29;
    % read flood asci file
    filename_i = sprintf('%s%d.asc',filename, year_start+e_i-1);
    fprintf('%s\n',filename_i);
    ascifile = [climada_global.data_dir filesep 'hazards' filesep 'Flood_Barisal' filesep 'today' filesep filename_i];
    flood_grid = flipud(dlmread(ascifile,delimiter,R1,C1));
    nodata_value = -9999;
    flood_grid(flood_grid==nodata_value) = 0;

    % write hazard intensity into structure
    hazard.intensity(e_i ,:) = reshape(flood_grid,1,no_col*no_row);
    hazard.name{1,e_i}       = filename_i;
    hazard.yyyy(1,e_i)       = year_start+e_i-1;
    hazard.datenum(1,e_i)    = datenum(hazard.yyyy(e_i),1,1);
end

%save in climada/data/hazard
hazard_filename = [climada_global.data_dir filesep 'hazards' filesep 'BCC_hazard_FL.mat'];
save(hazard_filename, 'hazard')
fprintf('saved hazard in %s\n',hazard_filename);


%%





    



