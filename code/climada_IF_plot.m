function climada_IF_plot(hazard,logscale_check,schematic_check)
% visualize Annual Expected Intensity per centroid as a map
% NAME:
%   climada_IF_plot
% PURPOSE:
%   plot annual expected intensities
% CALLING SEQUENCE:
%   climada_IF_plot(hazard, logscale_check, schematic_check)
% EXAMPLE:
%   climada_IF_plot(hazard, 1, 1)
% INPUTS:
%   hazard:     standard climada hazard set structure
% OPTIONAL INPUT PARAMETERS:
%   schematic_check:    whether to plot axis labels, colorbar (default = 0)
%                       or not (= 1)
%   logscale_check:     whether to use a logarithmic colorbar
% OUTPUTS:
%   figure
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150703 init
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('hazard'                  ,'var'),    hazard          =[];    end
if ~exist('logscale_check'          ,'var'),    logscale_check  = 1;   	end
if ~exist('schematic_check'         ,'var'),    schematic_check = 0;   	end

% PARAMETERS
% prompt for event damage set if not given
if isempty(hazard) % local GUI
    hazard=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard for intensity-frequency visualisation:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard=fullfile(pathname,filename);
    end
end

for centroid_i = 1:length(hazard.centroid_ID)
    hazard.IF_at_centroid(centroid_i,1) = full(sum(hazard.intensity(:,centroid_i) .* hazard.frequency'));
end

% load the event damage set, if a filename has been passed
if ~isstruct(hazard)
    hazard_file=hazard;hazard=[];
    load(hazard_file);
end

% create the figure
scale  = max(hazard.lon) - min(hazard.lon);
ax_buffer = 3; %ax_buffer = 30;
ax_lim = [min(hazard.lon)-scale/ax_buffer  max(hazard.lon)+scale/ax_buffer ...
    max(min(hazard.lat),-60)-scale/ax_buffer  min(max(hazard.lat),80)+scale/ax_buffer];

markersize = 2;

% fig = climada_figuresize(height,height*scale2+0.15);
cmap = climada_colormap('schematic');
nz = hazard.IF_at_centroid>0;
% cbar = plotclr(hazard.lon(nz), hazard.lat(nz), hazard.IF_at_centroid(nz),'s',markersize,1,...
%     [],[],colormap(cmap),0,logscale_check);

[x, y] = meshgrid(unique(hazard.lon),unique(hazard.lat));
gridded_h_FI_at_c = griddata(hazard.lon,hazard.lat,hazard.IF_at_centroid,x,y);
gridded_h_FI_at_c(gridded_h_FI_at_c==0) = nan;
contourf(x,y,gridded_h_FI_at_c, 'edgecolor','none'); hold on
cbar = colorbar;
colormap(cmap);

if logscale_check
    caxis(log([min(hazard.IF_at_centroid(nz)) max(hazard.IF_at_centroid(nz))]))
else
    caxis([min(hazard.IF_at_centroid(nz)) max(hazard.IF_at_centroid(nz))])
end

if ~isfield(hazard,'reference_year'),   hazard.reference_year = climada_global.present_reference_year; end

name_str = sprintf('Annual expected intensities for %s',num2str(hazard.reference_year));

cbar_lbl = 'Expected Intensity (%s)';
if logscale_check
    label_str = sprintf([cbar_lbl ' (exponential scale)'],hazard.units);
else
    label_str = sprintf(cbar_lbl,hazard.units);
end
set(get(cbar,'ylabel'),'String', label_str ,'fontsize',12);

if schematic_check
    set(gca,'xtickLabel',[],'ytickLabel',[])
    set(cbar,'Location','East','Visible','off')
else
    xlabel('Longitude')
    ylabel('Latitude')
end

box on
climada_plot_world_borders(1)
axis(ax_lim)
axis equal
axis(ax_lim)
title({name_str,strrep(hazard.comment,'_',' ')})