function [data, x, y] = climada_grid2array(data_grid, reference_box)

data = zeros(numel(data_grid),1);
tmp_x = zeros(numel(data_grid),1);
tmp_ = zeros(numel(data_grid),1);

[size_y, size_x] = size(data_grid);

if size_x <= size_y
    tmp_grid = data_grid;
    trans_check = 0;
else
    tmp_grid = data_grid';
    trans_check = 1;
end

[size_y, size_x] = size(tmp_grid);
dx = (reference_box(2) - reference_box(1))/size_x;
dy = (reference_box(4) - reference_box(3))/size_y;

t0 = clock;
for i = 1 : size_x
    ndx = (i-1)*size_y;
    data(ndx + 1 : ndx + size_y)    = tmp_grid(i,:);
    tmp_y(ndx + 1 : ndx + size_y,1) = (size_x - i+1) * dy;
    tmp_x(ndx + 1 : ndx + size_y,1) = (1:size_y) .* dx;
end
time_elapsed = clock - t0;
sprintf('converting from grid to array took %f sec',time_elapsed);

if trans_check
    x = tmp_y; y = tmp_x;
else
    x = tmp_x; y = tmp_y;
end

x = x + reference_box(1);
y = y + reference_box(3);