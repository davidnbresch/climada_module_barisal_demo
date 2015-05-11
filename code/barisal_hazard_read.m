% barisal_hazard_read

hazard_names = {'flood_depth_monsoon' 'flood_depth_cyclone' 'flood_duration_monsoon' 'flood_duration_cyclone' 'cyclone_wind'}; 

%climate change scenario
cc_scenario = {'no change' 'moderate' 'extreme'}; 

%timehorizon
timehorizon = [2014 2030 2050];

% set foldername
foldername = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\risk_modelling\flood\';


%% loop over hazards
for h_i = 1:length(hazard_names)
    hazard_name = hazard_names{h_i};
    
    switch hazard_name
        case 'flood_depth_monsoon';
            filename = 'MonsoonMaxInundationDepths1983.asc';
            units = 'm';
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_depth_monsoon'];
            
        case 'flood_depth_cyclone';
            filename = 'CycloneMaxInundationDepths1988.asc';  
            units = 'm';
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_depth_cyclone'];
            
        case 'flood_duration_monsoon';
            filename = 'MonsoonMaxInundationDurations1983.asc';
            units = 'days';
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration_monsoon'];
            
        case 'flood_duration_cyclone';
            filename = 'CycloneMaxInundationDurations1988.asc';
            units = 'days';
            hazard_set_file = [climada_global.data_dir filesep 'hazards' filesep 'Barisal_BCC_hazard_FL_duration_cyclone'];
    end   
    
    % loop over climate change scenarios
    for cc_i = 1:length(cc_scenario)
        switch cc_scenario{cc_i}
            case 'no change'
                folder_ = 'Floods_current_scenario\';
                reference_year = 2014;
                hazard_set_file_ = sprintf('%s_%d.mat',hazard_set_file,2014);
                
            case 'moderate'
                folder_ = 'Floods_2030_CCMod\';
                reference_year = 2030;
                hazard_set_file_ = sprintf('%s_cc_%d_%s.mat',hazard_set_file,reference_year,cc_scenario{cc_i});
                 
            case 'extreme'
                folder_ = 'Floods_2030_CCHigh\';
                reference_year = 2030;
                hazard_set_file_ = sprintf('%s_cc_%d_%s.mat',hazard_set_file,reference_year,cc_scenario{cc_i});
        end
        
        hazard = climada_asci2hazard([foldername folder_ filename]);
        hazard.comment = sprintf('Modelled by W+B, %s, %s', hazard_name, cc_scenario{cc_i});
        hazard.reference_year = reference_year;
        hazard.peril_ID = 'FL';
        hazard.units    = units;
        if strcmp(hazard_name,'flood_depth_cyclone')|strcmp(hazard_name,'flood_duration_cyclone')
            hazard.frequency = [0.05 0.05 0.02 0.03]*1.2;
        end
        
        save(hazard_set_file_, 'hazard');        
    end

end

%%





