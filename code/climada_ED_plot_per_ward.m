function fig = climada_ED_plot_per_ward(EDS,entity, BCC_wards, timehorizon, hazard_name)
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

% assign values to wards
ED_per_ward = zeros(length(BCC_wards),1);
for w_i = 1:length(BCC_wards)
    indx = entity.assets.Ward == w_i;
    ED_per_ward(w_i) = sum(EDS(t_i).ED_at_centroid(indx));
end
fprintf('\t - Total damage for %s in %d is BDT %2.0f''000\n',hazard_name,timehorizon, sum(ED_per_ward))


fig = climada_figuresize(0.75,0.75);
no_colors    = 10;
cbar         = jet(no_colors);
values       = ED_per_ward;
pos_indx     = values>0;
min_value    = min(values(pos_indx));
max_value    = max(values);
range_values = linspace(min_value,max_value,no_colors);

BCC_ward_no  = [BCC_wards.Ward_no];
for ward_i = 1:length(BCC_ward_no) 
    indx   = find(values(ward_i)<=range_values);
    indx   = indx(1);
    indx_w = find(BCC_ward_no == ward_i);
    h      = fill(BCC_wards(indx_w).lon,BCC_wards(indx_w).lat,cbar(indx,:));
    hold on
    ward_label = sprintf('Ward %d', BCC_wards(indx_w).WARDS_F_ID);
    text(entity.assets.lon(ward_i), entity.assets.lat(ward_i), ward_label,...
            'Horizontalalignment','center','verticalalignment','bottom','color','k') %brighten(cbar(indx,:),-0.8)
    %text(entity.assets.lon(ward_i), entity.assets.lat(ward_i), sprintf('Ward %d', ward_i),...
    %        'Horizontalalignment','center','verticalalignment','bottom','color','k') %brighten(cbar(indx,:),-0.8)
    g = plot(entity.assets.lon(ward_i), entity.assets.lat(ward_i),'kx','markersize',8,'LineWidth',1.5);    
end
colormap(cbar)
t = colorbar;
%cbar_label = sprintf('Intensity %s (%s)', hazard.peril_ID, hazard.units);
set(get(t,'ylabel'),'String', ('1000 BDT'),'fontsize',12);
caxis([min_value max_value])
%axislim = [min(EDS(1).assets.lon) max(EDS(1).assets.lon)*1 min(EDS(1).assets.lat) max(EDS(1).assets.lat)*1];
%axislim = [min(hazard.lon) max(hazard.lon)*1 min(hazard.lat) max(hazard.lat)*1];
%axislim = [90.25 90.45 22.6 22.8]; %barisal close up BCC 
axislim = [90.297 90.3957 22.64 22.752]; %barisal close up BCC 
axis(axislim)
axis equal
%titlestr = sprintf('%d, Annual damage, %s - %s', timehorizon(t_i), EDS(1).annotation_name, EDS(end).annotation_name);
titlestr = sprintf('%d: Annual damage from %s', timehorizon, strrep(hazard_name,'_',' '));
title({titlestr})
legend(g,'Lat/lon coordinates for assets')


