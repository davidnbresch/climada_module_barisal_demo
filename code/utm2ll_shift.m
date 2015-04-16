function [lon, lat] = utm2ll_shift(X,Y,datum)
% climada
% MODULE:
%   barisal_demo
% NAME:
%   utm2ll_shift
% PURPOSE:
%   Transform local UTM/btm coordinates (GCS  Everest 1830, used in
%   Barisal, Bangladesh, X,Y (in meters)) into worldwide
%   lon-lat-coordinates (in degrees). Additionally a shift is added to
%   improve the transformation. Default datum is WGS84.
% CALLING SEQUENCE:
%   [lon,lat] = utm2ll_shift(X,Y,datum)
% EXAMPLE:
%   [lon,lat] = utm2ll_shift(X,Y)
% INPUTS:
%   X: Eastin meter
%   Y: North in meter
% OPTIONAL INPUT PARAMETERS:
%   DATUM can be astring in the following list:
%		'wgs84': World Geodetic System 1984 (default)
%		'nad27': North American Datum 1927
%		'clk66': Clarke 1866
%		'nad83': North American Datum 1983
%		'grs80': Geodetic Reference System 1980
%		'int24': International 1924 / Hayford 1909
% OUTPUTS:
% MODIFICATION HISTORY:
% Lea Mueller, muellele@gmail.com, 20150315, initial
% Lea Mueller, muellele@gmail.com, 20150415, order change, first lon, then lat
%-

lon=[]; lat=[]; % init output

% poor man's version to check arguments
% and to set default value where  appropriate
if ~exist('X','var'),return;end
if ~exist('Y','var'),return;end
if ~exist('datum','var'), datum = ''; end


% PARAMETERS
% manual shift to correct offset (~270m easting, 305m northing)
delta_X = -270;
delta_Y = +305;

% Calculation
[lon,lat] = utm2ll(X+delta_X, Y+delta_Y, datum);





    
    
    

