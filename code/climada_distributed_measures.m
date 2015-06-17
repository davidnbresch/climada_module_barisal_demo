function [hazard,measure]= climada_distributed_measures(measure_file, hazard,bsxfun_op,save_file,check_plot)
%climada_distributed_measures
% NAME:
%   climada_distributed_measures
% PURPOSE:
%  function to account for measures that are introduced via an ascii file
%  the script reduces a selected hazard intensity by the grid value after
%  checking for NaN value and ensuring no minus values are produced.
% CALLING SEQUENCE:
%   climada_distributed_measures
% EXAMPLE:
%   climada_distributed_measures
% INPUTS:
%   acsii grid
%   .mat hazard file
% OPTIONAL INPUT PARAMETERS:
% OUTPUTS:
%   new hazard file called hazard_w_measures, hazard, measure_grid (reduction from ascii
%   file)
% RESTRICTIONS:
%   none
% MODIFICATION HISTORY:
% Jacob Anz, j.anz@gmx.net, 20150609
% Gilles Stassen, gillesstassen@hotmail.com, 20150610, clean up, bsx_fun_op
%                   added
%-

global climada_global
if ~climada_init_vars,return;end % init/import global variables

if ~exist('measure_file',       'var'), measure_file = '';      end
if ~exist('hazard',             'var'), hazard = [];            end
if ~exist('check_plot',         'var'), check_plot = '';        end
if ~exist('save_file',          'var'), save_file = '';         end
if ~exist('bsxfun_op',          'var'), bsxfun_op = [];         end

%open the grid measure (ascii) file
[fP, fN, fE] = fileparts(measure_file);
[~,fN_]      = fileparts(fP);

if strcmp(fE,'.asc') || isempty(fE)
    measure=climada_ascii_read(measure_file);
elseif strcmp(fE,'.mat')
    load(measure_file)
end
if ~exist('measure','var')
    cprintf([1 0 0],'ERROR: invalid measure file\n'); return
elseif ~isstruct(measure)
    cprintf([1 0 0],'ERROR: invalid measure file\n'); return
elseif ~isfield(measure,'lon') || ~isfield(measure,'lat') || ~isfield(measure,'value')
    cprintf([1 0 0],'ERROR: invalid measure file\n'); return
end
    
%change the ccordinates from UTM to degrees
%check if coordinates are already in UTM
if numel(num2str(abs(measure.lon(1))))>8
    [measure.lon,measure.lat]=utm2ll_shift(measure.lon,measure.lat);
end

if isempty(hazard) % local GUI
    hazard=[climada_global.data_dir filesep 'hazards' filesep '*.mat'];
    [filename, pathname] = uigetfile(hazard, 'Select hazard event set (or centroids) to encode to:');
    if isequal(filename,0) || isequal(pathname,0)
        return; % cancel
    else
        hazard=fullfile(pathname,filename);
    end
end
hazard_ori = hazard;

% load the hazard, if a filename has been passed
if ~isstruct(hazard_ori)
    hazard_file=hazard_ori;hazard_ori=[];
    load(hazard_file);
end

%Check if the the reduction contains 'ans', if it does, turn it to 0 or in
%case of a fraction turn it to 1

F_measure = scatteredInterpolant(measure.lon',measure.lat',measure.value');

measure.value_at_centroids = F_measure(hazard_ori.lon,hazard_ori.lat);

if isempty(bsxfun_op)
    if min(measure.value_at_centroids) >= -1 && max(measure.value_at_centroids) <= 2
        % probably multiplication factor
        bsxfun_op = @times;
    elseif any(measure.value_at_centroids <0) || (max(abs(measure.value_at_centroids)) >2)
        % probably addition
        bsxfun_op = @plus;
    else
        %Ask if the file is a fraction or a total value
        choice = menu('Select if the ascii file contains total values (e.g. "m/s") or fractions ("%")','Total value','Fraction');
        if choice==1
            bsxfun_op = @plus;
        else
            bsxfun_op = @times;
        end
    end
end

% set the NaNs to identity for specified operator
switch func2str(bsxfun_op)
    case 'plus'
        measure.value_at_centroids(isnan(measure.value_at_centroids))=0;
    case 'times'
        measure.value_at_centroids(isnan(measure.value_at_centroids))=1;
    case 'divide'
        measure.value_at_centroids(isnan(measure.value_at_centroids))=1;
end

% calculate impact on hazard
hazard.intensity = bsxfun(bsxfun_op,hazard_ori.intensity,measure.value_at_centroids);

% if choice ==1;
%     measure.value_2(isnan(measure.value_2))=0;
%     for event_i=1:length(hazard.event_ID)
%         hazard_w_measures.intensity(event_i,:)=full(hazard.intensity(event_i,:))+ measure.value_at_centroids;
%     end  
% else
%     measure.value_2(isnan(measure.value_2))=1;
%     for event_i=1:length(hazard.event_ID)
%         for centroid_i=1:length(hazard.centroid_ID)
%             hazard_w_measures.intensity(event_i,centroid_i)=full(hazard.intensity(event_i,centroid_i)).* measure.value_at_centroids(centroid_i);
%         end
%     end
% end


%check for negative intensities and set them to 0
hazard.intensity(hazard.intensity <0) = 0;
% [find1,find2]=find(full(hazard_w_measures.intensity<0));
% hazard_w_measures.intensity(find1,find2)=0;

measure.comment = strrep([fN_],'_',' ');
hazard.comment= strrep([fN_ ' ' fN],'_',' ');
if ~isempty(save_file)
    hazard.filename = save_file;
end

%plot the difference of the hazard sets without and with measures

hazard_diff=hazard_ori;
hazard_diff.intensity=hazard.intensity-hazard_ori.intensity;
if check_plot
    figure('color','w')
    climada_hazard_plot_hr(hazard_diff);
    title(sprintf('Impact of %s on %s hazard %d',measure.comment,hazard.peril_ID,hazard.reference_year))
    figure('color','w'); hold on
    s = scatter(measure.lon,measure.lat,'filled');
    set(s,'cdata',measure.value,'marker','s')
    climada_plot_world_borders(1.5)
end