function [data, x, y] = climada_grid2array(data_grid, x_vector, y_vector)

data = [];      x = [];     y = [];

if ~exist('data_grid',      'var'),     return;             end
if ~exist('x_vector',       'var'),     x_vector = [];      end
if ~exist('y_vector',       'var'),     y_vector = [];      end

if numel(x_vector) == 4 && isempty(y_vector)
    reference_box = x_vector;
elseif ~isempty(y_vector) && ~isempty(x_vector)
    reference_box = [min(x_vector) max(x_vector) min(y_vector) max(y_vector)];
else
    reference_box = [];
end

switch ndims(data_grid)
    case 2
        data    = zeros(numel(data_grid),1);
        tmp_x   = zeros(numel(data_grid),1);
        tmp_y   = zeros(numel(data_grid),1);
        
        [size_y, size_x] = size(data_grid);
        
        if isempty(reference_box)
            reference_box = [1 size_x 1 size_y];
        end
        
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
        
        %t0 = clock;
        for i = 1 : size_x
            ndx = (i-1)*size_y;
            data    (ndx + 1 : ndx + size_y)    = tmp_grid(:,i);
            tmp_y   (ndx + 1 : ndx + size_y,1)  = i*dy;%(size_x - i+1) * dy;
            tmp_x   (ndx + 1 : ndx + size_y,1)  = (1:size_y) .* dx;
        end
        %time_elapsed = clock - t0;
        %fprintf('converting from grid to array took %f sec \n',time_elapsed);
        
        if trans_check
            x = tmp_y; y = tmp_x;
        else
            x = tmp_x; y = tmp_y;
        end
        
        x = x + reference_box(1);
        y = y + reference_box(3);
    case 3
        [size_y, size_x, size_t] = size(data_grid);
        
        try
            data    = zeros(size_x*size_y,size_t);
        catch
            data    = sparse(size_x*size_y,size_t);
        end
        tmp_x   = zeros(size_x*size_y,1);
        tmp_y   = zeros(size_x*size_y,1);
        
        if isempty(reference_box)
            reference_box = [1 size_x 1 size_y];
        end
        
        if size_x <= size_y
            tmp_grid = data_grid;
            trans_check = 0;
        else
            tmp_grid = permute(data_grid,[2 1 3]);
            trans_check = 1;
        end
        
        [size_y, size_x,size_t] = size(tmp_grid);
        dx = (reference_box(2) - reference_box(1))/size_x;
        dy = (reference_box(4) - reference_box(3))/size_y;
        
        for i = 1 : size_x
            ndx = (i-1)*size_y;
            data    (ndx + 1 : ndx + size_x,:)  = squeeze(tmp_grid(i,:,:));
            tmp_y   (ndx + 1 : ndx + size_y,1)  = (size_x - i+1) * dy;
            tmp_x   (ndx + 1 : ndx + size_y,1)  = (1:size_y) .* dx;
        end
        
        if trans_check
            x = tmp_y; y = tmp_x;
        else
            x = tmp_x; y = tmp_y;
        end
        
        x = x + reference_box(1);
        y = y + reference_box(3);
    otherwise
        return;
end