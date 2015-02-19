function [precip_grid, precip_array] = climada_import_global_precip(precip_nc_file,centroids_rect,time_frame,check_plot)
% http://www.esrl.noaa.gov/psd/data/gridded/data.gpcc.html

precip_grid = []; precip_array = [];

if ~exist('precip_nc_file','var') || ~exist(precip_nc_file, 'file')
    [fN, fP] = uigetfile('*.nc', 'Select precipitation data file');
    if isequal(fN,0) || isequal(fP,0),  return;             end
    precip_nc_file = [fP fN];
end

if ~exist('centroids_rect','var'),      centroids_rect = [];    end
if ~exist('time_frame','var'),          time_frame = [];        end
if ~exist('check_plot','var'),          check_plot = 0;         end

fprintf('reading precipitation data...')
precip_nc_info = ncinfo(precip_nc_file);
for fld_i =  1 : numel(precip_nc_info.Variables)
    fld_name = precip_nc_info.Variables(fld_i).Name;
    precip_grid.(fld_name)=ncread(precip_nc_file,fld_name);
end
fprintf('done \n')

if ~all(ismember({'lon' 'lat' 'time' 'precip'},fieldnames(precip_grid)))
    cprintf([1 0.5 0],'WARNING: data missing, try again with a different .nc file \n');
    return;
end

if any(precip_grid.lon > 180)
    precip_grid.lon(precip_grid.lon > 180) = precip_grid.lon(precip_grid.lon > 180)- 360;
    precip_grid.precip = [precip_grid.precip(precip_grid.lon < 0,:,:); precip_grid.precip(precip_grid.lon >= 0,:,:)];
    precip_grid.lon = [precip_grid.lon(precip_grid.lon < 0,:,:); precip_grid.lon(precip_grid.lon >= 0,:,:)];
end

precip_grid.info = precip_nc_info.Attributes;
precip_grid.file = precip_nc_info.Filename;

precip_grid.start_day = datestr(datenum(1901,1,1) + precip_grid.time(1),'dd/mm/yyyy');
precip_grid.time = precip_grid.time + datenum(1901,1,1);

if ~isempty(centroids_rect)
    lon_logical = ceil(precip_grid.lon) >= centroids_rect(1) & floor(precip_grid.lon) <= centroids_rect(2);
    lat_logical = ceil(precip_grid.lat) >= centroids_rect(3) & floor(precip_grid.lat) <= centroids_rect(4);
    precip_grid.lon = precip_grid.lon(lon_logical);
    precip_grid.lat = precip_grid.lat(lat_logical);
    precip_grid.precip = precip_grid.precip(lon_logical,lat_logical,:);
end

if ~isempty(time_frame)
    time_logical = precip_grid.time >= time_frame(1) & precip_grid.time <= time_frame(2);
    if any(time_logical)
        precip_grid.time = precip_grid.time(time_logical);
        precip_grid.precip = precip_grid.precip(:,:,time_logical);
    else
        fprintf('WARNING: invalid time frame, reverting to full data series')
    end
end

precip_grid.precip(precip_grid.precip<0) = 0;

if check_plot
    % plot rainfall for month with heaviest rainfall in time frame
    [~,max_rain_month] = max(sum(sum(precip_grid.precip,2),1));
    figure('name', 'Most intense rainfall');
    title( sprintf('Precipitation for the month of %s in mm',datestr(precip_grid.time(max_rain_month),'mm/yyyy')));
    hold on
    imagesc(precip_grid.lon,precip_grid.lat,permute(precip_grid.precip(:,:,max_rain_month),[2 1 3]));
    set(gca,'YDir','normal');
    climada_plot_world_borders(2)
    axis equal
    axis (centroids_rect);
    c_map = [linspace(1,0,12)' linspace(1,0,12)' linspace(1,0.75,12)'];
    colormap(c_map)
    colorbar
    box on
    grid on
    hold off
    
    figure('name',sprintf('Average monthly precipitation for period %s to %s',datestr(precip_grid.time(1)), datestr(precip_grid.time(end))),'Color','w');
    title(sprintf('Average monthly precipitation for period %s to %s',datestr(precip_grid.time(1)), datestr(precip_grid.time(end))));
    hold on
    for month_i = 1 : 12
        try
            monthly_mean_rainfall(month_i) = ...
                mean(mean(mean(precip_grid.precip(:,:,month(precip_grid.time) == month_i),1),2));
            monthly_stdv_rainfall(month_i) = ...
                std(mean(mean(precip_grid.precip(:,:,month(precip_grid.time) == month_i),1),2));
        catch
            monthly_mean_rainfall(month_i) = ...
                mean(squeeze(mean(mean(precip_grid.precip(:,:,str2double(datestr((precip_grid.time),'mm')) == month_i),1),2)));
            monthly_stdv_rainfall(month_i) = ...
                std(squeeze(mean(mean(precip_grid.precip(:,:,str2double(datestr((precip_grid.time),'mm')) == month_i),1),2)));
        end
    end
    bar(1:12,monthly_mean_rainfall, 'facecolor', [0.0   0.4078   0.5451])
    hold on
    errorbar(1:12,monthly_mean_rainfall,monthly_stdv_rainfall,'.k')
    colormap(c_map)
    set(gca,'Ylim',[0 ceil((max(monthly_mean_rainfall)+max(monthly_stdv_rainfall))/100)*100]);
    set(gca,'Xtick',1:12)
    set(gca,'XTickLabel', {'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'});
    xlabel('Month')
    ylabel('Rainfall [mm]')
    hold off
end

if nargout == 2
    [precip_array.precip, precip_array.lon, precip_array.lat] = ...
        climada_grid2array(precip_grid.precip, [min(precip_grid.lon) max(precip_grid.lon) min(precip_grid.lat) max(precip_grid.lat)]);
    precip_array.time       = precip_grid.time;
    precip_array.info       = precip_grid.info;
    precip_array.start_day  = precip_grid.start_day;
    precip_array.file       = precip_grid.file;
end


