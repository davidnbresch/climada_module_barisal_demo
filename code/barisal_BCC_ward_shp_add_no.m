
%% assign ward numbers to ward shapefile

% load BCC wards
BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards.mat'];
% BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards_Ward_no_added.mat'];
load(BCC_wards_savename)
indx2 = strfind(BCC_wards_savename,filesep);
fprintf('\t - loaded BCC specifics: %s\n', BCC_wards_savename(indx2(end)+1:end))

% load entity
[hazard, entity, label] = barisal_hazard_entity_load('flood_depth', 'no change', 2010);

ward_no = unique(entity.assets.Ward_Nr);
for w_i=1:length(ward_no)
    indx = find(entity.assets.Ward_Nr == ward_no(w_i));
    ward_points(w_i,1) = mean(entity.assets.lon(indx));
    ward_points(w_i,2) = mean(entity.assets.lat(indx));
end
% ward_points = [entity.assets.lon(1:30) entity.assets.lat(1:30)];

% go through all ward points and find polygon which is point is inside
for w_i=1:length(BCC_wards)
    polygon_nodes = [BCC_wards(w_i).lon' BCC_wards(w_i).lat'];
    [cn,on] = inpoly(ward_points,polygon_nodes);
    cn = find(cn);
    if ~isempty(cn)
        BCC_wards(w_i).Ward_no = cn(1);
    else
        BCC_wards(w_i).Ward_no = 0;
    end
end

% hand corrections
% BCC_wards(2).Ward_no = 6;
% BCC_wards(3).Ward_no = 5;
% BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards.mat'];
% save(BCC_wards_savename,'BCC_wards')

BCC_wards( 2).Ward_no = 31;
BCC_wards(13).Ward_no =  7;
BCC_wards( 3).Ward_no = 34;
BCC_wards( 6).Ward_no = 32;
BCC_wards_savename = [climada_global.data_dir filesep 'entities' filesep 'BCC_wards_no_added.mat'];
save(BCC_wards_savename,'BCC_wards')


% check with figure
climada_figuresize(0.7,0.5)
for w_i=1:length(BCC_wards)
    hold on
    h(3)= plot(BCC_wards(w_i).lon,BCC_wards(w_i).lat,'color',[244 164 96 ]/255);%sandybrown;
    %text(mean(BCC_wards(w_i).lon), mean(BCC_wards(w_i).lat), int2str(BCC_wards(w_i).WARDS_F_ID))
    text(mean(BCC_wards(w_i).lon), mean(BCC_wards(w_i).lat), sprintf('%d:%d',w_i, BCC_wards(w_i).Ward_no))
end
plot(ward_points(:,1),ward_points(:,2), 'x','markersize',2)
for a_i = 1:length(ward_no)
    text(ward_points(a_i,1),ward_points(a_i,2), int2str(ward_no(a_i)),'color','b','Horizontalalignment','center','verticalalignment','bottom')
end





