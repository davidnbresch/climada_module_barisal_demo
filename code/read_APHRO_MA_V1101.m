function [precip_grid, precip_array] = read_APHRO_MA_V1101(years, centroids_rect, precip_save_file, location, check_plot)
% climada TC hazard event set generate
% NAME:
%   climada_tc_hazard_set
% PURPOSE:
%   Read in the Monsoon Asia daily precipitation data from the APHRODITE
%   gridded .nc data file. The original data was downloaded from 
%   http://www.chikyu.ac.jp/precip/cgi-bin/aphrodite/script/aphrodite_cgi.cgi/download?file=%2FV1101%2FAPHRO_MA_V1101%2F025deg%2Fnc
%   The information is structered into both a gridded structure, and a
%   singleton array structure, which may then be further processed.
%   next: climada_ma_hazard_set
%   
%   See the following for more information on the original source:
%   http://www.chikyu.ac.jp/precip/cgi-bin/aphrodite/script/aphrodite_cgi.cgi/download?file=%2FV1101%2FAPHRO_MA_V1101%2F025deg%2Fnc
%   http://iridl.ldeo.columbia.edu/SOURCES/.RIHN/.aphrodite/.V1003R1/
%   http://www.chikyu.ac.jp/precip/cgi-bin/aphrodite/script/aphrodite_cgi.cgi/register
%
% CALLING SEQUENCE:
%   [precip_grid, precip_array] = read_APHRO_MA_V1101(years, centroids_rect, precip_save_file, location, check_plot)
% EXAMPLE:
%   [precip_grid, precip_array] = read_APHRO_MA_V1101(1995:2001, centroids_rect, [], location, 1)
%   [~, precip_array] = read_APHRO_MA_V1101(1995:2001, centroids_rect, [], location, 1)
% INPUTS:
%   years:          a vector containing years of interest
%   centroids_rect: row vector of the form [min_lon max_lon min_lat max_lat] 
%                   defining the area of interest, to be cropped out of
%                   original data set.
% OPTIONAL INPUT PARAMETERS:
%   precip_save_file:   path and filename of save location for output files
%   location:           struct with fields .longitude, .latitude, which
%                       defines a 'sampling point' for graph of rainfall timeseries
%   check_plots:        whether to show time lapse of daily data;
% OUTPUTS:
%   precip_grid:        the data in gridded form. Struct with fields
%                       .precip: NxMxT array with precip amounts [mm] along
%                           two spatial dimensions, M, N, and time dimension, T.
%                       .time:   1xT vector with temporal coords
%                       .lon:    1xN vector with longitudinal coords of grid
%                       .lon:    1xM vector with latitudinal coords of grid
%   precip_array:       the data in singleton array form. Struct with fields
%                       .precip: (NxM)xT array with precip amounts [mm] along
%                           one spatial dimension, NxM, and time dimension T.
%                       .time:   1xT vector with temporal coords
%                       .lon:    NxM vector with longitudinal coords for each point
%                       .lon:    NxM vector with latitudinal coords for each point
% MODIFICATION HISTORY:
% Gilles Stassen, gillesstassen@hotmail.com, 20150121
%-

precip_grid.time = [];      precip_array.time = [];
precip_grid.lon = [];       precip_array.lon = [];
precip_grid.lat = [];       precip_array.lat = [];
precip_grid.precip = [];    precip_array.precip = [];

if ~climada_init_vars, return; end
global climada_global;

if ~exist('years',              'var'), years = [];             end
if ~exist('centroids_rect',     'var'), centroids_rect = [];    end
if ~exist('location',           'var'), location = [];          end
if ~exist('check_plot',         'var'), check_plot = 0;         end
if ~exist('precip_save_file',   'var'), precip_save_file = [];  end

if isempty(years)
    prompt   ='Choose years of interest (separate using commas):';
    name     ='Years ';
    default_ans = {'1993'};
    answer = inputdlg(prompt,name,1,default_ans);
    answer = cell2mat(answer);
    years = sscanf(answer,'%i,');
    
    if max(years) > 2007 || min(years) < 1951
        years_ignored = [years(years<1951); years(years > 2007)];
        years = years(~ismember(years, years_ignored));
        years_ignored = num2str(years_ignored');
        cprintf([1 0.5 0],'WARNING: data only available from 1951 to 2007, ignoring input years %s \n',years_ignored);
    end
end
years = sort(years, 'ascend');
root_dir = fileparts(pwd);

format_str = '%s';

for year_i = 1 : numel(years)
    year = num2str(years(year_i));
    
    APHRO_nc_file = ['APHRO_MA_025deg_V1101.' year '.nc'];
    
    try
        APHRO_nc_file = subdir([APHRO_dir filesep APHRO_nc_file]);
    catch
        APHRO_nc_file = subdir([root_dir filesep APHRO_nc_file]);
        APHRO_dir = fileparts(APHRO_nc_file.name);
    end
    
    if isempty(APHRO_nc_file)
        cprintf([1 0.5 0],'WARNING: APHRO data file not found for specified year \n')
        [fN, fP] = uigetfile('*.nc', 'Select precipitation data file');
        if isequal(fN,0) || isequal(fP,0),  return;             end
        APHRO_nc_file = [fP fN];
    end
    
    msg_str = sprintf('reading daily precipitation data for %s... ',year);
    fprintf(format_str,msg_str);
    format_str=[repmat('\b',1,length(msg_str)) '%s'];
    
    precip_nc_info = ncinfo(APHRO_nc_file.name);
    for fld_i =  1 : numel(precip_nc_info.Variables)
        fld_name = precip_nc_info.Variables(fld_i).Name;
        tmp_p_g.(fld_name)=ncread(APHRO_nc_file.name,fld_name);
    end

    
    % crop to centroids rect
    if ~isempty(centroids_rect)
        lon_logical = ceil(tmp_p_g.longitude) >= centroids_rect(1) & floor(tmp_p_g.longitude) <= centroids_rect(2);
        lat_logical = ceil(tmp_p_g.latitude) >= centroids_rect(3) & floor(tmp_p_g.latitude) <= centroids_rect(4);
        tmp_p_g.longitude = tmp_p_g.longitude(lon_logical);
        tmp_p_g.latitude = tmp_p_g.latitude(lat_logical);
        tmp_p_g.precip = tmp_p_g.precip(lon_logical,lat_logical,:);
    end
    
    tmp_p_g.precip(tmp_p_g.precip < 0) = 0;
    tmp_p_g.precip = double(tmp_p_g.precip);
    
    [tmp_p_a.precip, tmp_p_a.longitude, tmp_p_a.latitude] = ...
        climada_grid2array(tmp_p_g.precip, tmp_p_g.longitude, tmp_p_g.longitude);
    tmp_p_a.time = tmp_p_g.time;
    
    precip_grid.time(end+1:end+numel(tmp_p_g.time)) = tmp_p_g.time + datenum(year, 'yyyy');
    precip_grid.precip(:,:,end+1:end+numel(tmp_p_g.time)) = tmp_p_g.precip;    
    precip_array.precip(:,end+1:end+numel(tmp_p_a.time))  = tmp_p_a.precip;
end

fprintf('done \n')

precip_grid.precip(:,:,1) = [];

precip_array.time = precip_grid.time;
precip_grid.lon = tmp_p_g.longitude;
precip_grid.lat = tmp_p_g.latitude;
precip_array.lon = tmp_p_a.longitude;
precip_array.lat = tmp_p_a.latitude;

% precip_grid.time = sort(precip_grid.time,'ascend');
% precip_grid.precip = sort(precip_grid.precip,3,'ascend');
% precip_array.time = sort(precip_array.time,'ascend');
% precip_array.precip = sort(precip_array.precip,2,'ascend');

if check_plot
    n_days = numel(tmp_p_g.time);
    figure
    hold on
    set(gca,'YDir','normal');
    map_tr = [linspace(1,0,12)' linspace(1,0,12)' linspace(1,0.75,12)'];
    colormap(map_tr)
    climada_plot_world_borders;
    axis equal
    axis(centroids_rect)
    %axis tight
    view(0,90)
    hold on
    for day_i = 1 : n_days
        title(sprintf('Daily precipitation in %s [mm]',datestr(datenum(year,'yyyy')+day_i)));
        if any(any(tmp_p_g.precip(:,:,day_i)))
            s = surf(tmp_p_g.lon,tmp_p_g.lat,permute(tmp_p_g.precip(:,:,day_i)-100,[2 1 3]));
            shading('interp')
            pause(0.2)
            delete(s);
        else
            pause(0.2)
        end
        % contourf(precip_grid.longitude,precip_grid.latitude,permute(precip_grid.precip(:,:,day_i),[2 1 3]));
        % imagesc(precip_grid.longitude,precip_grid.latitude,permute(precip_grid.precip(:,:,day_i),[2 1 3]));
    end
    hold off
end

if isstruct(location)
    figure('name', sprintf('Daily precipitation in %s during %s',strrep(location.name,' ',''),year), 'color', 'w');
    hold on
    title(sprintf('Daily precipitation in %s during %s',strrep(location.name,' ',''),year));
    r = climada_geo_distance(location.longitude, location.latitude,precip_array.lon, precip_array.lat);
    [~, r_ndx] = min(r);
    plot(precip_array.time, precip_array.precip(r_ndx,:))
    ylabel('Precipitation [mm]')
    xlabel('Date')
    set(gca,'xlim',[min(precip_array.time) max(precip_array.time)])
    set(gca,'xticklabel',datestr(get(gca,'xtick'),'mm/yyyy'))
    hold off
end

if ~isempty(precip_save_file) && ischar(precip_save_file)
    save(precip_save_file, 'precip_grid', 'precip_array')
end

return

