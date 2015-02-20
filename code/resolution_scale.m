function [x_s, y_s] = resolution_scale(x,y,x_factor,y_factor)
% climada resolution_upscale
% NAME:
%   resolution_upscale
% PURPOSE:
%   Scale a set of x, y axes according to scaling factors. Specifically,
%   the function is designed to increase the resolution of a grid specified
%   by vectors x, y of the form:
%           x = [1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4]
%           y = [1 2 3 4 1 2 3 4 1 2 3 4 1 2 3 4]
%   (arranged vice versa also works)
% CALLING SEQUENCE:
%   [x_hr, y_hr] = resolution_upscale(x,y,x_factor,y_factor)
% EXAMPLE:
%   [x_hr, y_hr] = resolution_upscale(x,y,x_factor,y_factor)
% INPUTS:
%   x               : original x data
%   y               : original y data
% OPTIONAL INPUT PARAMETERS:
%   x_factor        : scaling factor for the x direction
%   y_factor        : scaling factor for the y direction
% OUTPUTS:
%   x_s             : scaled x data
%   y_s             : scaled y data
% MODIFICATION HISTORY:
%   Gilles Stassen, gillesstassen@hotmail.com 20141104
%   Gilles Stassen, gillesstassen@hotmail.com 20150108 - improved speed ~x3
%-
x_s = x; y_s = y;

err_str = 'WARNING: One or more non-optional inputs not provided - data not scaled';

if ~exist('x','var'),           fprintf(err_str);   return;     end
if ~exist('y','var'),           fprintf(err_str);   return;     end
if ~exist('x_factor','var') &&...
        ~exist('y_factor','var'),fprintf(err_str);  return;     end
if ~exist('y_factor','var'),    y_factor = 1;       return;     end
if ~exist('x_factor','var'),    x_factor = 1;       return;     end

x_res = (max(x) - min(x))/(numel(unique(x))-1);
y_res = (max(y) - min(y))/(numel(unique(y))-1);

x_s_res = x_res/x_factor;
y_s_res = y_res/y_factor;

x_s_size = (numel(unique(x))-1)*x_factor + 1;
y_s_size = (numel(unique(y))-1)*y_factor + 1;

x_s = zeros(1,x_s_size); y_s = zeros(1,y_s_size);

for i = 1 : x_s_size
    ndx = (i-1)*y_s_size;
    y_s(1,ndx + 1 : ndx + y_s_size) = (x_s_size - i+1) * y_s_res    + min(y);
    x_s(1,ndx + 1 : ndx + y_s_size) = (1:y_s_size) .* x_s_res       + min(x);
end

