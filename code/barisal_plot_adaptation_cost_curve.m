
% plot selection of measures impact for barisal
% Lea Mueller, 20160531, for ECA publication Barisal

measures_filename = [climada_global.data_dir filesep 'results' filesep 'measures_impact_barisal.xlsx'];
[measures, measures_impact] = climada_measures_read(measures_filename);
measures_impact.color_keep = 1;



fig = climada_figuresize(0.4,0.95);
climada_adaptation_cost_curve(measures_impact)
print(fig,'-dpdf',[climada_global.data_dir filesep 'results' filesep 'adaptation_cost_curve_barisal.pdf'])


% climada_adaptation_cost_curve(measures_impact,'','','','','','',1)