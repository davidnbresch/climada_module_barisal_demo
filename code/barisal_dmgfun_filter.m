function [DamFunID indx_depth indx_duration indx_windspeed] = barisal_dmgfun_filter(entity)

% find damagefunctions for flood duration (max intensity value is 14 days)
% fprintf('\t- Organize damage functions (flood depth, flood duration, cyclone wind speed)\n')
max_value_duration =  14.0;
max_value_depth    =   4.0;
max_value_windspeed= 400.0;

DamFunID       = unique(entity.damagefunctions.DamageFunID);
indx_depth     = zeros(size(DamFunID));
indx_duration  = zeros(size(DamFunID));
indx_windspeed = zeros(size(DamFunID));
for d_i = 1:length(DamFunID)
    indx = entity.damagefunctions.DamageFunID == DamFunID(d_i);
    max_value = max(entity.damagefunctions.Intensity(indx));
    
    if max_value<=max_value_depth+10^-2 & max_value>=max_value_depth-10^-2
        indx_depth(d_i) = 1;
    elseif max_value<=max_value_duration+10^-2 & max_value>=max_value_duration-10^-2
        indx_duration(d_i) = 1;
    elseif max_value<=max_value_windspeed+10^-2 & max_value>=max_value_windspeed-10^-2
        indx_windspeed(d_i) = 1;
    end
end
indx_depth     = logical(indx_depth);
indx_duration  = logical(indx_duration);
indx_windspeed = logical(indx_windspeed);
