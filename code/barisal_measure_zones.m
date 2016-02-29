



filename = 'M:\BGCC\CHR\RK\RS\A_Sustainable_Development\Projects\ECA\BarisalBangladesh\risk_modelling\4_measures\adaptation_strategy\Measure_zone.shp';
zones = shaperead(filename);

% figure
% climada_shapeplotter(zones,'Zones_code')

figure
fill(zones(1).X(~isnan(zones(1).X)),zones(1).Y(~isnan(zones(1).X)),'-b')
hold on
fill(zones(2).X(~isnan(zones(2).X)),zones(2).Y(~isnan(zones(2).X)),'-c')
fill(zones(3).X(~isnan(zones(3).X)),zones(3).Y(~isnan(zones(3).X)),'-g')
fill(zones(4).X(~isnan(zones(4).X)),zones(4).Y(~isnan(zones(4).X)),'-r')
for s = 1:4
    text(mean(zones(s).X(~isnan(zones(s).X))),mean(zones(s).Y(~isnan(zones(s).X))),int2str(s),'fontsize',14)
end
axis equal



resilient_buildings_zone_B.X = zones(1).X(~isnan(zones(1).X));
resilient_buildings_zone_B.Y = zones(1).Y(~isnan(zones(1).Y));
[resilient_buildings_zone_B.lon resilient_buildings_zone_B.lat] = utm2ll_shift(resilient_buildings_zone_B.X, resilient_buildings_zone_B.Y);


figure
fill(resilient_buildings_zone_B.X, resilient_buildings_zone_B.Y ,'-b')

figure
fill(resilient_buildings_zone_B.lon, resilient_buildings_zone_B.lat ,'-b')
