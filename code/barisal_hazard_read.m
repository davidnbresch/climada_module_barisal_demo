% barisal_hazard_read
module_data_dir=[climada_modules_dir filesep 'climada_module_barisal_demo' filesep 'data'];
spec = '_rain_only';

hazard_names = {'flood_depth_monsoon' 'flood_depth_cyclone' 'flood_duration_monsoon' 'flood_duration_cyclone'};% 'cyclone_wind'}; 

%climate change scenario
cc_scenario = {'no change' 'moderate' 'extreme'}; 

%timehorizon
timehorizon = [2014 2030 2050];
future_reference_year = 2050;

% set foldername
% foldername = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\risk_modelling\flood\';
% foldername = [climada_global.data_dir filesep 'hazards' filesep 'Flood_Barisal\'];
foldername = [module_data_dir filesep 'hazards' filesep 'Flood_Barisal' filesep];
foldername = '\\CHRB1048.CORP.GWPNET.COM\homes\E\S1T2EN\Documents\climada_GIT\vulnerability_analysis_barisal\3_results\';

%% loop over hazards
for y_i = future_reference_year
for h_i = 1:length(hazard_names)
    hazard_name = hazard_names{h_i};
    
    switch hazard_name
        case 'flood_depth_monsoon';
            filename = 'MonsoonMaxInundationDepths1983.asc';
            units = 'm';
            hazard_set_file = [module_data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_depth_monsoon' ];
            
        case 'flood_depth_cyclone';
            filename = 'CycloneMaxInundationDepths10-11-2007.asc';  
            units = 'm';
            hazard_set_file = [module_data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_depth_cyclone' ];
            
        case 'flood_duration_monsoon';
            filename = 'MonsoonMaxInundationDurations1983.asc';
            units = 'days';
            hazard_set_file = [module_data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration_monsoon' ];
            
        case 'flood_duration_cyclone';
            filename = 'CycloneMaxInundationDurations10-11-2007.asc';
            units = 'days';
            hazard_set_file = [module_data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration_cyclone' ];
    
    end   
    
    % loop over climate change scenarios
    for cc_i = 1:length(cc_scenario)
        switch cc_scenario{cc_i}
            case 'no change'
                folder_ = ['Floods_current_scenario' filesep];
                reference_year = 2014;
                hazard_set_file_ = sprintf('%s_%d.mat',hazard_set_file,reference_year);
                
            case 'moderate'
                folder_ = ['Floods_' num2str(y_i) '_CCMod' filesep];
                reference_year = y_i;
                %folder_ = 'Floods_2050_CCMod\';
                %reference_year = 2050;
                hazard_set_file_ = sprintf('%s_cc_%d_%s.mat',hazard_set_file,reference_year,cc_scenario{cc_i});
                 
            case 'extreme'
                folder_ = ['Floods_' num2str(y_i) '_CCHigh' filesep];
                reference_year = y_i;
                %folder_ = 'Floods_2050_CCHigh\';
                %reference_year = 2050;
                hazard_set_file_ = sprintf('%s_cc_%d_%s.mat',hazard_set_file,reference_year,cc_scenario{cc_i});
        end
        
        fprintf('%s\n',[foldername folder_ filename])
%         hazard = climada_asci2hazard([foldername folder_ filename],6);
        hazard = climada_asci2hazard([foldername folder_ filename]);

        hazard.comment = sprintf('%s %s, modelled by W+B', strrep(hazard_name,'_',' '), cc_scenario{cc_i});
        hazard.reference_year = reference_year;
        if ~isempty(strfind(upper(hazard_name), 'DURATION'))
            hazard.peril_ID = 'FL_';
        else
            hazard.peril_ID = 'FL';
        end
        
        hazard.units    = units;
        if strfind(hazard_name,'cyclone')
            hazard.frequency    = [0.05 0.05 0.02 0.03]*1.2;
            hazard.dd           = [10 22 25 30];
            hazard.mm           = [11 11 05 11];
            hazard.yyyy         = [2007 1998 2009 1988];
            hazard.datenum      = datenum(hazard.yyyy,hazard.mm,hazard.dd);
        end
        
        save(hazard_set_file_, 'hazard');        
    end
end
end

%%





