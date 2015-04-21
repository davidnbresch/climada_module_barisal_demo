function screw = barisal_climate_screw(peril_ID, time_horizon, scen_type)
% NAME:
% barisal_climate_screw
% PURPOSE:
%   Generate a climate change screw structure for the city of Barisal,
%   based on the data and methodology outlined in the ECA Barisal base line
%   report.
% CALLING SEQUENCE:
%   screw = barisal_climate_screw(peril_ID, time_horizon, scen_type)
% EXAMPLE:
%   screw = barisal_climate_screw('TC', 2050, 'extreme')
% INPUTS:
%   if no input provided, the default climate screw for TC 2030 moderate
%   scenario will be returned.
% OPTIONAL INPUT PARAMETERS:
%   peril_ID:       defines the peril to which the screw applies
%   time_horizon:   defines the year at which the screw applies [2030 or 2050]
%   scen_type:      defines the scenario type to which the screw applies [moderate or extreme]
% OUTPUTS:
%   screw:      defines the climate change scenario. A 1xN structure with
%               fields:
%                   .hazard_fld     defines the hazard field to be changed
%                   .change         extent of the change at time horizon
%                   .year           time horizon
%                   .hazard_crit    hazard field to which criteria apply
%                   .criteria       criteria for events/locations to change
%                   .bsxfun_op      operation of change (e.g. @times,@plus) (function handle)
%               specifying N transformations to the original hazard set.
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com, 20150421 
%-

if ~exist('peril_ID'    ,'var'),    peril_ID        = 'TC';         end
if ~exist('time_horizon','var'),    time_horizon    = 2030;         end
if ~exist('scen_type'   ,'var') || ~ismember(scen_type, {'moderate' 'extreme'})
    scen_type       = 'moderate';   
end
       

% screw = 4x4 struct: [ 2030 moderate   2050 moderate
%                       2030 extreme    2050 extreme   ]

screw_IDs = {   '2030_moderate'   '2050_moderate'
                '2030_extreme'    '2050_extreme'   };
screw_ndx = strcmp(screw_IDs, sprintf('%i_%s',time_horizon,scen_type));

%% TC  
% TC 2030 moderate
barisal_screw(1,1).TC(1).hazard_fld = 'intensity';
barisal_screw(1,1).TC(2).hazard_fld = 'frequency';
barisal_screw(1,1).TC(1).change     =  1.006;
barisal_screw(1,1).TC(2).change     =  0.991;
barisal_screw(1,1).TC(1).year       =  2030;
barisal_screw(1,1).TC(2).year       =  2030;
barisal_screw(1,1).TC(1).hazard_crit= 'category';
barisal_screw(1,1).TC(2).hazard_crit= 'category';
barisal_screw(1,1).TC(1).criteria   =  [3 4 5];
barisal_screw(1,1).TC(2).criteria   =  [1 2 3 4 5];
barisal_screw(1,1).TC(1).bsxfun_op  =  @times;
barisal_screw(1,1).TC(2).bsxfun_op  =  @times;

% TC 2050 moderate
barisal_screw(1,2).TC(1).hazard_fld = 'intensity';
barisal_screw(1,2).TC(2).hazard_fld = 'frequency';
barisal_screw(1,2).TC(1).change     =  1.010;
barisal_screw(1,2).TC(2).change     =  0.982;
barisal_screw(1,2).TC(1).year       =  2050;
barisal_screw(1,2).TC(2).year       =  2050;
barisal_screw(1,2).TC(1).hazard_crit= 'category';
barisal_screw(1,2).TC(2).hazard_crit= 'category';
barisal_screw(1,2).TC(1).criteria   =  [3 4 5];
barisal_screw(1,2).TC(2).criteria   =  [1 2 3 4 5];
barisal_screw(1,2).TC(1).bsxfun_op  =  @times;
barisal_screw(1,2).TC(2).bsxfun_op  =  @times;

% TC 2030 extreme
barisal_screw(2,1).TC(1).hazard_fld = 'intensity';
barisal_screw(2,1).TC(2).hazard_fld = 'frequency';
barisal_screw(2,1).TC(1).change     =  1.040;
barisal_screw(2,1).TC(2).change     =  1.070;
barisal_screw(2,1).TC(1).year       =  2030;
barisal_screw(2,1).TC(2).year       =  2030;
barisal_screw(2,1).TC(1).hazard_crit= 'category';
barisal_screw(2,1).TC(2).hazard_crit= 'category';
barisal_screw(2,1).TC(1).criteria   =  [3 4 5];
barisal_screw(2,1).TC(2).criteria   =  [1 2 3 4 5];
barisal_screw(2,1).TC(1).bsxfun_op  =  @times;
barisal_screw(2,1).TC(2).bsxfun_op  =  @times;

% TC 2050 extreme
barisal_screw(2,2).TC(1).hazard_fld = 'intensity';
barisal_screw(2,2).TC(2).hazard_fld = 'frequency';
barisal_screw(2,2).TC(1).change     =  1.060;
barisal_screw(2,2).TC(2).change     =  1.130;
barisal_screw(2,2).TC(1).year       =  2050;
barisal_screw(2,2).TC(2).year       =  2050;
barisal_screw(2,2).TC(1).hazard_crit= 'category';
barisal_screw(2,2).TC(2).hazard_crit= 'category';
barisal_screw(2,2).TC(1).criteria   =  [3 4 5];
barisal_screw(2,2).TC(2).criteria   =  [1 2 3 4 5];
barisal_screw(2,2).TC(1).bsxfun_op  =  @times;
barisal_screw(2,2).TC(2).bsxfun_op  =  @times;

%% TS (subsidence + sea level rise)
% TS 2030 moderate
barisal_screw(1,1).FL.hazard_fld    = 'intensity';
barisal_screw(1,1).FL.change        =  0.14;
barisal_screw(1,1).FL.year          =  2030;
barisal_screw(1,1).FL.bsxfun_op     =  @plus;

% TS 2050 moderate
barisal_screw(1,2).FL.hazard_fld    = 'intensity';
barisal_screw(1,2).FL.change        =  0.33;
barisal_screw(1,2).FL.year          =  2050;
barisal_screw(1,2).FL.bsxfun_op     =  @plus;

% TS 2030 extreme
barisal_screw(2,1).FL.hazard_fld    = 'intensity';
barisal_screw(2,1).FL.change        =  0.30;
barisal_screw(2,1).FL.year          =  2030;
barisal_screw(2,1).FL.bsxfun_op     =  @plus;

% TS 2050 extreme
barisal_screw(2,2).FL.hazard_fld    = 'intensity';
barisal_screw(2,2).FL.change        =  0.53;
barisal_screw(2,2).FL.year          =  2050;
barisal_screw(2,2).FL.bsxfun_op     =  @plus;

%% RF (monsoon rainfall intensity)
% RF 2030 moderate
barisal_screw(1,1).RF.hazard_fld    = 'intensity';
barisal_screw(1,1).RF.change        =  1.020;
barisal_screw(1,1).RF.year          =  2030;
barisal_screw(1,1).RF.hazard_crit   = 'mm';
barisal_screw(1,1).RF.criteria      =  [6 7 8 9];
barisal_screw(1,1).RF.bsxfun_op     =  @times;

% RF 2050 moderate
barisal_screw(1,2).RF.hazard_fld    = 'intensity';
barisal_screw(1,2).RF.change        =  0.040;
barisal_screw(1,2).RF.year          =  2050;
barisal_screw(1,2).RF.hazard_crit   = 'mm';
barisal_screw(1,2).RF.criteria      =  [6 7 8 9];
barisal_screw(1,2).RF.bsxfun_op     =  @times;

% RF 2030 extreme
barisal_screw(2,1).RF.hazard_fld    = 'intensity';
barisal_screw(2,1).RF.change        =  0.050;
barisal_screw(2,1).RF.year          =  2030;
barisal_screw(2,1).RF.hazard_crit   = 'mm';
barisal_screw(2,1).RF.criteria      =  [6 7 8 9];
barisal_screw(2,1).RF.bsxfun_op     =  @times;

% RF 2050 extreme
barisal_screw(2,2).RF.hazard_fld    = 'intensity';
barisal_screw(2,2).RF.change        =  0.070;
barisal_screw(2,2).RF.year          =  2050;
barisal_screw(2,2).RF.hazard_crit   = 'mm';
barisal_screw(2,2).RF.criteria      =  [6 7 8 9];
barisal_screw(2,2).RF.bsxfun_op     =  @times;

%assign correct screw to output
screw = barisal_screw(screw_ndx).(peril_ID);

end