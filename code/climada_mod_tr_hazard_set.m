function hazard = climada_mod_tr_hazard_set(hazard_tr, hazard_precip,location, hazard_set_file, check_plot)
% climada
% NAME:
%   climada_mod_tr_hazard_set
% PURPOSE:
%   a validation routine which compares the rainfields generated
% CALLING SEQUENCE:
%   hazard=climada_mod_tr_hazard_set(hazard_tr, hazard_ma,location, hazard_set_file, check_plot)
% EXAMPLE:
%   hazard=climada_mod_tr_hazard_set(hazard_tr, hazard_ma,[], hazard_set_file, 1)
% INPUTS:
%   hazard_tr:          original TR hazard set as generated using the
%                       R-CLIPER rainfield. See climada_tr_hazard_set
%   hazard_precip:      reference precipitation data for same region as
%                       hazard_tr, with climada hazard structure
%   location:           can either be a struct with fields .longitude and
%                       .latitude, in which case the nearest centroid is
%                       found to specified coords, or a Nx1 array with
%                       centroid_IDs, defining 'test' points
%   hazard_set_file: the file path and name of the modified TR hazard set.
% OPTIONAL INPUT PARAMETERS:
%   check_plots:        whether to show IFC plot comparison
% OUTPUTS:
%   hazard: a hazard event set, see core climada doc
%       also written to a .mat file (see hazard_set_file)
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmsil.com, 20150121
% TO DO:    generalise the modification routine - numerical approach to
%           finding best itnensity and frequency modifiers.
%-

hazard =[]; % init

global climada_global
if ~climada_init_vars,return;end % init/import global variables

% poor man's version to check arguments
if ~exist('hazard_tr',          'var'), hazard_tr           =[];    end
if ~exist('hazard_set_file',    'var'), hazard_set_file     ='';    end
if ~exist('hazard_precip',      'var'), hazard_precip       =[];    end
if ~exist('check_plot',         'var'), check_plot           =0;    end
if ~exist('location',           'var'),
    if isfield(hazard_tr, 'centroid_ID')
        centroid_ID = median(hazard_tr.centroid_ID);
    elseif isfield(hazard_precip, 'centroid_ID')
        centroid_ID = median(hazard_precip.centroid_ID);
    else
        fprintf('ERROR: no field .centroid_ID in hazard struct \n')
        return;
    end
else
    if isstruct(location)
        dist2loc = climada_geo_distance(location.lon,location.lat,hazard_tr.lon,hazard_tr.lat);
        [~,min_dist_ndx] = min(dist2loc);
        centroid_ID = hazard_tr.centroid_ID(min_dist_ndx);
    elseif isvector(location)
        centroid_ID = location;
    else
        if isfield(hazard_tr, 'centroid_ID')
            cprintf([1 0.5 0],'WARNING: invalid input location, reverting to default centroid \n');
            centroid_ID = median(hazard_tr.centroid_ID);
        elseif isfield(hazard_precip, 'centroid_ID')
            cprintf([1 0.5 0],'WARNING: invalid input location, reverting to default centroid \n');
            centroid_ID = median(hazard_precip.centroid_ID);
        else
            fprintf('ERROR: can not find valid centroid_ID \n')
            return;
        end
    end
end

% hazard args passed in opposite order
if strcmp(hazard_precip.peril_ID,'TR')
    hazard = hazard_precip;
    hazard_precip = hazard_tr;
    hazard_tr = hazard;
else
    hazard = hazard_tr;
end
hazard.peril_ID = 'TR_m';

freq_mod = [linspace(15,1/1000000,numel(hazard.frequency))];
%freq_mod = [linspace(1,1/1000000,numel(hazard.frequency))];
[~, sort_ndx] = sort(hazard.intensity(:,centroid_ID),'ascend');
freq_mod = freq_mod(sort_ndx);

%int_mod = linspace(1/4,1/10,numel(hazard.frequency));
int_mod = linspace(1/2,1/6,numel(hazard.frequency));
%int_mod = ones(size(hazard.frequency));
int_mod = int_mod(sort_ndx);
int_mod = repmat(int_mod',1,numel(hazard.centroid_ID));

int_offset = 0;
freq_offset = 0;

hazard.frequency = (hazard.frequency.*freq_mod)+freq_offset;
hazard.intensity = (int_mod.*hazard.intensity)+int_offset;

if check_plot
    IFC_precip = climada_hazard2IFC(hazard_precip,centroid_ID);
    IFC_tr = climada_hazard2IFC(hazard_tr,centroid_ID);
    IFC = climada_hazard2IFC(hazard,centroid_ID);
    figure('color', 'w')
    hold on
    title('Intensity-frequency curves for MA, TR and modified TR hazard sets')
    climada_IFC_plot(IFC_precip,1,1,1,1,1);
    hold on
    climada_IFC_plot(IFC,1,1,1,1,2);
    hold on
    climada_IFC_plot(IFC_tr,1,1,1,1,3);
    hold off
end

if ~isempty(hazard_set_file)
    fprintf('saving modified TR hazard set to %s \n', hazard_set_file)
    save(hazard_set_file,'hazard')
end


return
