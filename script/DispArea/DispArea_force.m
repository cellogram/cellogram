function [area, TF, in] = DispArea(V, disp, z_thres, silent)
    force = ones (205,3);
    if (nargin<3)
        z_thres = -1;
    end
    if (nargin<4)
        silent = false;
    end
    gap = 5;
    xx = round(min(V(:, 1)))-gap:round(max(V(:, 1)))+gap;
    yy = round(min(V(:, 2)))-gap:round(max(V(:, 2)))+gap;
    [x, y] = meshgrid(xx, yy);
    % try the following meshgrid, in order to have um for sure / check if
    % it gives the same area
    %[x,y] = meshgrid (1:1:max(cellogram.mesh.points(:,1))*scaling, 1:1:max(cellogram.mesh.points(:,2))*scaling)
    query = [x(:), y(:)];
    RBF_res = RBFinterp(V(:, 1:2), disp, query, 'GS', 3);
    
    % find the query location with largest z-deformation
    t = RBF_res(:, 3);
    idx = find(t == min(t));
    mx = floor(idx / length(yy)) + 1;
    my = idx - (mx-1) * length(yy);
    mx = xx(mx);
    my = yy(my);
    
    % BFS to grow area
    delta_xy = [0, 1; 1, 0; 0, -1; -1, 0];
    query_in_area = [mx, my];
    to_search = [mx, my];
    visited = [ToIdx(mx, my)];

    while size(to_search, 1) > 0
        for dir = 1:4
            new_xy = to_search(1, :) + delta_xy(dir, :);
            if new_xy(1)<xx(1), continue;end
            if new_xy(1)>xx(end), continue;end
            if new_xy(2)<yy(1), continue;end
            if new_xy(2)>yy(end), continue;end
            if ismember(ToIdx(new_xy(1), new_xy(2)), visited)
                continue;
            end
            visited = [visited; ToIdx(new_xy(1), new_xy(2))];
            idx = length(yy) * (new_xy(1)-xx(1)) + (new_xy(2)-yy(1)) + 1;
            if RBF_res(idx, 3) > z_thres
                continue;
            end
            
            % new loc
            query_in_area = [query_in_area; new_xy];  %#ok
            to_search = [to_search; new_xy];  %#ok
        end
        
        to_search = to_search(2:end, :);
    end
    
    % add fake force vector
    force = rand(205,3);
    scaling = 1; % should be scaling = cellogram.analysis_settings.scaling;
    k = boundary(query_in_area(:,1),query_in_area(:,2),1);
    p = polyshape(query_in_area(k,1),query_in_area(k,2));
    % result
    if (~silent)
        figure
        hold on
        quiver3(query(:, 1), query(:, 2), zeros(length(RBF_res), 1), RBF_res(:, 1), RBF_res(:, 2), RBF_res(:, 3))
        scatter3(query_in_area(:, 1), query_in_area(:, 2), zeros(size(query_in_area, 1), 1));
        quiver3(V(:, 1), V(:, 2), V(:, 3), disp(:, 1), disp(:, 2), disp(:, 3))
        str = string(z_thres); title(str);
        hold on
        plot(p)
    end
    area = size(query_in_area, 1);
    
    in = inpolygon(query(:, 1),query(:, 2),p.Vertices(:,1),p.Vertices(:,2)); % which coordinates of the full grid are inside the area
    tmag = sqrt(force(:,1).^2 + force(:,2).^2 + force(:,3).^2 ); %xyz
    [TF] = griddata(V(:,1),V(:,2),tmag,query(:, 1),query(:, 2)); % interpolate the force values to the generated full grid
    mean(TF); % gives NaN, not all coordinates get a force value
    mean(TF(in)) % this should be correct - gives mean traction force over the affected area
    
end


function [res] = ToIdx(x, y)
    res = x*10000 + y;
end
