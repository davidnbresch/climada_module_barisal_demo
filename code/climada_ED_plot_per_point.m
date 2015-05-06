function fig = climada_ED_plot_per_point(EDS,entity, BCC_wards, timehorizon, hazard_name)
% create figure to plot expected damage per ward in barisal
% MODULE: 
%   barisal_demo
% NAME:
%   climada_ED_plot_per_ward
% PURPOSE:
%   plot annual expected damage per ward on a map
% CALLING SEQUENCE:
%   fig = climada_ED_plot_per_ward(EDS,entity, BCC_wards, timehorizon, hazard_name)
% EXAMPLE:
%   fig = climada_ED_plot_per_ward(EDS,entity, BCC_wards, timehorizon, hazard_name)
% INPUTS:
%   EDS: event damage set, as e.g. returned by climada_EDS_calc or
%       a file containing such a structure
%       EDS can contain multiple EDS, however default is that first EDS
%       will be used.
%   entity: an entity (see climada_entity_read)
%       > promted for if not given
%   BCC_wards: structure with shape file (polygons for all 30 wards in
%       Barisal)
% OPTIONAL INPUT PARAMETERS:
%   timehorizon: just for figure title, empty if not specified
%   hazard_name: just for figure title, empty if not specified
% OUTPUTS:
%   figure with damage per ward on a map 
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150429, init
%-


global climada_global
if ~climada_init_vars, return; end

% poor man's version to check arguments
if ~exist('EDS'           ,'var'), EDS         = []; end
if ~exist('entity'        ,'var'), entity      = []; end
if ~exist('BCC_wards'     ,'var'), BCC_wards   = []; end
if ~exist('timehorizon'   ,'var'), timehorizon = ''; end
if ~exist('hazard_name'   ,'var'), hazard_name = ''; end

fig = []; % init

if isempty(EDS),return;end 
if isempty(entity),return;end
if isempty(BCC_wards),return;end


% Set time horizon 
% If EDS contains multiple entries, t_i defines the EDS to be selected. 
% Default set to 1.
t_i = 1;


% find unique lat lons
[lon_lat,indx, indx2] = unique([EDS.assets.lon EDS.assets.lat],'rows');
lon_lat_all           = [EDS.assets.lon EDS.assets.lat];

% assign values to wards
ED_per_point = zeros(length(lon_lat),1);
for p_i = 1:length(ED_per_point)
    indx3 = find(indx2 == p_i);
    ED_per_point(p_i) = sum(EDS(t_i).ED_at_centroid(indx3));
end
fprintf('\t - Total damage at %d lon/lat-points %s in %d is BDT %2.0f''000\n',length(lon_lat),hazard_name,timehorizon, sum(ED_per_point))

% cmap = climada_colormap('damage');
fig = climada_figuresize(0.75,0.75);
no_colors    = 10;
cmap         = jet(no_colors);
values       = ED_per_point;
pos_indx     = values>0;
min_value    = min(values(pos_indx));
max_value    = max(values);
% range_values = linspace(min_value,max_value,no_colors);
mav = max_value*1.0;
markersize = 5;
% [h h_points] = plotclr(x,y,v, marker, markersize, colorbar_on, miv, mav, map, zero_off, v_exp)
[cbar,asset_handle]= plotclr(lon_lat(:,1), lon_lat(:,2), ED_per_point, 's',markersize, 1,min_value,max_value,cmap,0,0);
hold on
axis equal
xlabel('Longitude')
ylabel('Latitude')
set(gca,'layer','top')
set(get(cbar,'ylabel'),'string','Damage per point','fontsize',12)
climada_plot_world_borders(0.7);
axis equal
axislim = [90.285 90.3957 22.64 22.752]; %barisal close up BCC 
axis(axislim)
%titlestr = sprintf('%d, Annual damage, %s - %s', timehorizon(t_i), EDS(1).annotation_name, EDS(end).annotation_name);
titlestr = sprintf('%d: Annual damage from %s: BDT %2.0f''000\n', timehorizon, strrep(hazard_name,'_',' '), sum(ED_per_point));
title({titlestr})
box on
% legend(g,'Lat/lon coordinates for assets')


% loop over all words to plot color according to flood damage
BCC_ward_no  = [BCC_wards.Ward_no];
for ward_i = 1:length(BCC_ward_no) 
    %indx   = find(values(ward_i)<=range_values);
    %indx   = indx(1);
    %indx_w = find(BCC_ward_no == ward_i);
    plot(BCC_wards(ward_i).lon,BCC_wards(ward_i).lat,'color',[ 205 193 197 ]/255,'linewidth',2); %dark red
    text(mean(BCC_wards(ward_i).lon),mean(BCC_wards(ward_i).lat),BCC_wards(ward_i).UNION_NAME,...
        'Horizontalalignment','center','verticalalignment','bottom'); %grey
end








% % set the axis limits
% box on
% d = 0.015; %°Degree
% x_range = [min(entity.assets.lon)-d max(entity.assets.lon)+d];
% y_range = [min(entity.assets.lat)-d max(entity.assets.lat)+d];
% set(gca,'xlim',x_range,'ylim',y_range)

% colormap(cbar)
% t = colorbar;
% %cbar_label = sprintf('Intensity %s (%s)', hazard.peril_ID, hazard.units);
% set(get(t,'ylabel'),'String', ('1000 BDT'),'fontsize',12);
% caxis([min_value max_value])
% %axislim = [min(EDS(1).assets.lon) max(EDS(1).assets.lon)*1 min(EDS(1).assets.lat) max(EDS(1).assets.lat)*1];
% %axislim = [min(hazard.lon) max(hazard.lon)*1 min(hazard.lat) max(hazard.lat)*1];
% %axislim = [90.25 90.45 22.6 22.8]; %barisal close up BCC 


% % loop over all words to plot color according to flood damage
% BCC_ward_no  = [BCC_wards.Ward_no];
% for ward_i = 1:length(BCC_ward_no) 
%     indx   = find(values(ward_i)<=range_values);
%     indx   = indx(1);
%     indx_w = find(BCC_ward_no == ward_i);
%     h      = fill(BCC_wards(indx_w).lon,BCC_wards(indx_w).lat,cbar(indx,:));
%     hold on 
% end
% 
% % loop over all words to label ward name
% for ward_i = 1:length(BCC_ward_no) 
%     indx   = find(values(ward_i)<=range_values);
%     indx   = indx(1);
%     indx_w = find(BCC_ward_no == ward_i);
%     ward_label = BCC_wards(indx_w).UNION_NAME;
%     
%     text(lon_lat(ward_i,1), lon_lat(ward_i,2), ward_label,...
%             'Horizontalalignment','center','verticalalignment','bottom','color','k') %brighten(cbar(indx,:),-0.8)
%     %ward_label = sprintf('Ward %d', BCC_wards(indx_w).WARDS_F_ID);    
%     %text(entity.assets.lon(ward_i), entity.assets.lat(ward_i), ward_label,...
%     %        'Horizontalalignment','center','verticalalignment','bottom','color','k') %brighten(cbar(indx,:),-0.8)
%     %text(entity.assets.lon(ward_i), entity.assets.lat(ward_i), sprintf('Ward %d', ward_i),...
%     %        'Horizontalalignment','center','verticalalignment','bottom','color','k') %brighten(cbar(indx,:),-0.8)
%     %g = plot(entity.assets.lon(ward_i), entity.assets.lat(ward_i),'kx','markersize',8,'LineWidth',1.5);    
% end
% colormap(cbar)
% t = colorbar;
% %cbar_label = sprintf('Intensity %s (%s)', hazard.peril_ID, hazard.units);
% set(get(t,'ylabel'),'String', ('1000 BDT'),'fontsize',12);
% caxis([min_value max_value])
% %axislim = [min(EDS(1).assets.lon) max(EDS(1).assets.lon)*1 min(EDS(1).assets.lat) max(EDS(1).assets.lat)*1];
% %axislim = [min(hazard.lon) max(hazard.lon)*1 min(hazard.lat) max(hazard.lat)*1];
% %axislim = [90.25 90.45 22.6 22.8]; %barisal close up BCC 
% axislim = [90.285 90.3957 22.64 22.752]; %barisal close up BCC 
% axis(axislim)
% axis equal
% %titlestr = sprintf('%d, Annual damage, %s - %s', timehorizon(t_i), EDS(1).annotation_name, EDS(end).annotation_name);
% titlestr = sprintf('%d: Annual damage from %s: BDT %2.0f''000\n', timehorizon, strrep(hazard_name,'_',' '), sum(ED_per_ward));
% title({titlestr})
% % legend(g,'Lat/lon coordinates for assets')


