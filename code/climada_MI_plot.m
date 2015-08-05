function climada_MI_plot(EDS, percentage_of_value_flag,currency,unit_exp,logscale_check,schematic_check)
% visualize Annual Expected Damage per centroid as a map
% NAME:
%   climada_MI_plot
% PURPOSE:
%   plot annual expected measure impact
% CALLING SEQUENCE:
%   climada_MI_plot(EDS, percentage_of_value_flag)
% EXAMPLE:
%   climada_MI_plot(EDS, percentage_of_value_flag)
% INPUTS:
%   EDS output from climada_measures_impact_report (which has field
%   .MI_at_centroid)
% OPTIONAL INPUT PARAMETERS:
%   percentage_of_value_flag: Set to 1 if you wish to plot damages as
%                             percentage of asset values
% OUTPUTS:
%   figure
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150625 init
% Gilles Stassen, gillesstassen@hotmail.com, 20150628 - schematic_check added
% Lea Mueller, muellele@gmail.com, 20150706, add switch for UTM instead of lat/lon coordinates, inhibit limiting latitude values to -60° and +80°
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('EDS'                     ,'var'),    EDS     =   [];             end
if ~exist('percentage_of_value_flag','var'),    percentage_of_value_flag=0; end
if ~exist('currency'                ,'var'),    currency=   'USD';          end
if ~exist('unit_exp'                ,'var'),    unit_exp=   0;              end
if ~exist('logscale_check'          ,'var'),    logscale_check = 1;         end
if ~exist('schematic_check'         ,'var'),    schematic_check = 0;        end

% PARAMETERS
% prompt for event damage set if not given
if isempty(EDS) % local GUI
    EDS=[climada_global.data_dir filesep 'results' filesep '*.mat'];
    [filename, pathname] = uigetfile(EDS, 'Select EDS for ED visualisation:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        EDS=fullfile(pathname,filename);
    end
end

% load the event damage set, if a filename has been passed
if ~isstruct(EDS)
    EDS_file=EDS;EDS=[];
    load(EDS_file);
end

if ~isfield(EDS,'MI_at_centroid')
    cprintf([1 0 0],'ERROR: EDS must have field ''MI_at_centroid'', see climada_measures_impact_report\n')
    return
end

% set currency unit
if unit_exp == 0
    MI = sum([EDS.MI_at_centroid],1);
    while max(MI) > 1000
        unit_exp    = unit_exp + 3;
        MI          = MI./1000;
    end
end

unit_exp = interp1([0 3 6 9 12],[0 3 6 9 12],unit_exp,'nearest');
switch unit_exp
    case 0
        unit_char = '';
    case 3
        unit_char = 'k';
    case 6
        unit_char = 'm';
    case 9
        unit_char = 'bn';
    case 12
        unit_char = 'tn';
    otherwise
        unit_char = sprintf('10^%i',unit_exp);
end

for EDS_i = 1: length(EDS)
    if length(EDS) > 1
        figure
        hold on
    end
    if MI(EDS_i) == 0;
        cprintf([0 0 1],'NOTE: no benefit from measure %s, plot skipped\n',EDS(EDS_i).annotation_name)
        close; continue
    end
    % create the figure
    scale  = max(EDS(EDS_i).assets.lon) - min(EDS(EDS_i).assets.lon);
    ax_buffer = 10; %ax_buffer = 30;
    if max(EDS(EDS_i).assets.lon)>1000 %not lat/lon but UTM coordinates
        min_lat = 0;
        max_lat = 10^9;
    else %normal lat/lon coordinates, set maximum south to -60° and maximum north to 80°
        min_lat = -60;
        max_lat = 80;
    end
    ax_lim = [min(EDS(EDS_i).assets.lon)-scale/ax_buffer               max(EDS(EDS_i).assets.lon)+scale/ax_buffer ...
              max(min(EDS(EDS_i).assets.lat),min_lat)-scale/ax_buffer  min(max(EDS(EDS_i).assets.lat),max_lat)+scale/ax_buffer];
    
    markersize = 3;
    
    % to deal with multi-valued points
    lon_lat = unique([EDS(EDS_i).assets.lon EDS(EDS_i).assets.lat],'rows');
    if length(lon_lat) ~= length(EDS(EDS_i).assets.lon)
        for i = 1:length(lon_lat)
            ndx = EDS(EDS_i).assets.lon == lon_lat(i,1) & EDS(EDS_i).assets.lat == lon_lat(i,2);
            MI_sum_centroid(i)   = sum(EDS(EDS_i).MI_at_centroid(ndx));
            val_sum_centroid(i) = sum(EDS(EDS_i).assets.Value(ndx));
        end
    else
        lon_lat = [EDS(EDS_i).assets.lon EDS(EDS_i).assets.lat];
        MI_sum_centroid      = EDS(EDS_i).MI_at_centroid;
        val_sum_centroid    = EDS(EDS_i).assets.Value;
    end
    
    % fig = climada_figuresize(height,height*scale2+0.15);
    cmap_a = climada_colormap('benefit');
    cmap_b = flipud(climada_colormap('schematic'));
    if percentage_of_value_flag

        nz = MI_sum_centroid <=-1;
        dam_TAV = (MI_sum_centroid(nz) ./ val_sum_centroid(nz)) *100;
        plotclr(lon_lat(nz,1)', lon_lat(nz,2)', dam_TAV','s',markersize,0,...
            [],[],colormap(cmap_b),0,logscale_check);
        
        nz = MI_sum_centroid>=1 ;%| 
        dam_TAV = (MI_sum_centroid(nz) ./ val_sum_centroid(nz)) *100;
        cbar = plotclr(lon_lat(nz,1)', lon_lat(nz,2)', dam_TAV','s',markersize,1,...
            [],[],colormap(cmap_a),0,logscale_check);
        
        name_str = sprintf('Expected benefit (as percentage of value) for %s',num2str(EDS(EDS_i).reference_year));
    else
        nz = MI_sum_centroid <=-1;
        plotclr(lon_lat(nz,1), lon_lat(nz,2), MI_sum_centroid(nz),'s',markersize,0,...
            [],[],colormap(cmap_b),0,logscale_check);
        
        nz = MI_sum_centroid>=1 ;%| MI_sum_centroid <=-1;
        cbar = plotclr(lon_lat(nz,1), lon_lat(nz,2), MI_sum_centroid(nz),'s',markersize,1,...
            [],[],colormap(cmap_a),0,logscale_check);
        if logscale_check
            caxis(log([min(MI_sum_centroid(nz)) max(MI_sum_centroid(nz))]))
        else
            caxis([min(MI_sum_centroid(nz)) max(MI_sum_centroid(nz))])
        end
        name_str = sprintf('Expected benefit for %s:   %s %2.2f %s', ...
            num2str(EDS(EDS_i).reference_year),currency,...
           MI(EDS_i),unit_char);
        
    end
    if strfind(upper(currency),'PEOPLE')>0
        name_str = strrep(name_str,'benefit','no. of lives saved');
        name_str = strrep(name_str,'value','population');
        name_str = [strtok(name_str,':') ':' sprintf('   %2.2f %s %s',...
            MI(EDS_i),unit_char,currency)];
    end
    
    % set(fig,'Name',name_str)
    
    if percentage_of_value_flag
        cbar_lbl = 'Expected Benefit (% of value)';
        if strfind(upper(currency),'PEOPLE')>0,cbar_lbl = 'Expected no. of lives saved (% of population)'; end
        set(get(cbar,'ylabel'),'String',cbar_lbl,'fontsize',12);
    else
        cbar_lbl = 'Expected Benefit (%s)';
        if strfind(upper(currency),'PEOPLE')>0,cbar_lbl = 'Expected no. of lives saved'; end
        if logscale_check
            label_str = sprintf([cbar_lbl ' (exponential scale)'],currency);
        else
            label_str = sprintf(cbar_lbl,currency);
        end
        set(get(cbar,'ylabel'),'String', label_str ,'fontsize',12);
    end
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
    title({name_str,strrep(EDS(EDS_i).annotation_name,'_',' ')})
end
