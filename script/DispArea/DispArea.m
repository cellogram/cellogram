function [area] = DispArea(verts, disp, z_thres, silent)

    if (nargin<3)
        z_thres = -1;
    end
    if (nargin<4)
        silent = false;
    end
    gap = 5;
    xx = round(min(verts(:, 1)))-gap:round(max(verts(:, 1)))+gap;
    yy = round(min(verts(:, 2)))-gap:round(max(verts(:, 2)))+gap;
    [x, y] = meshgrid(xx, yy);
    query = [x(:), y(:)];
    RBF_res = RBFinterp(verts(:, 1:2), disp, query, 'GS', 3);
    
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
    
    % result
    if (~silent)
        figure
        hold on
        quiver3(query(:, 1), query(:, 2), zeros(length(RBF_res), 1), RBF_res(:, 1), RBF_res(:, 2), RBF_res(:, 3))
        scatter3(query_in_area(:, 1), query_in_area(:, 2), zeros(size(query_in_area, 1), 1));
        quiver3(verts(:, 1), verts(:, 2), verts(:, 3), disp(:, 1), disp(:, 2), disp(:, 3))
    end
    area = size(query_in_area, 1);
end


function [res] = ToIdx(x, y)
    res = x*10000 + y;
end
