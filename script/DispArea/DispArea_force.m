function [area, in, TF_mag_xyz, TF_mag_xy, TFx, TFy, TFz] = DispArea_force(V, disp, Tv, filepath, z_thres, silent)
%%   
% V = pts;
% z_thres = -3

     if (nargin<5)
         z_thres = -0.75;
     end
     if (nargin<6)
        silent = false;
     end
     %%
%     disp = D;
%     z_thres = -1;
    % for images that we start the z-stact from +X and finish at -Y !!!!!!!
    % check before use!!!
    disp(:, 3) = -disp(:, 3);
    
    gap = 1;
    xx = round(min(V(:, 1)))-gap:round(max(V(:, 1)))+gap;
    yy = round(min(V(:, 2)))-gap:round(max(V(:, 2)))+gap;
    [x, y] = meshgrid(xx, yy);
    % try the following meshgrid, in order to have um for sure / check if
    % it gives the same area
    %[x,y] = meshgrid (1:1:max(cellogram.mesh.points(:,1))*scaling, 1:1:max(cellogram.mesh.points(:,2))*scaling)
    query = [x(:), y(:)];
    RBF_res = RBFinterp(V(:, 1:2), disp, query, 'TPS2', 6);
    
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
    % force = rand(205,3);
    force = Tv;
    scaling = 1; % should be scaling = cellogram.analysis_settings.scaling;
    k = boundary(query_in_area(:,1),query_in_area(:,2),1);
    p = polyshape(query_in_area(k,1),query_in_area(k,2));
    % result
    if (~silent)
        fig = figure;
        hold on
        quiver3(query(:, 1), query(:, 2), zeros(length(RBF_res), 1), RBF_res(:, 1), RBF_res(:, 2), RBF_res(:, 3))
        %scatter3(query_in_area(:, 1), query_in_area(:, 2), zeros(size(query_in_area, 1), 1));
        plot(p)
        quiver3(V(:, 1), V(:, 2), V(:, 3), disp(:, 1), disp(:, 2), disp(:, 3))
        str_title = "z-thres = " + string(z_thres); title(str_title);
        legend("Interpolation", "Disp Area", "Raw Disp", "");
        grid on
        shrinkSize = 10;
        xlim([min(xx)+shrinkSize , max(xx)-shrinkSize ])
        ylim([min(yy)+shrinkSize , max(yy)-shrinkSize ])
    end
    
    filepath1 = filepath(1:end-1); FilePath = strrep(filepath, '\','/'); strrep(FilePath, '//','/'); slashIdx = strfind(FilePath, '/');
    pathSegment = FilePath(slashIdx(end-1)+1:slashIdx(end)-1);
    str = sprintf(pathSegment); title(str);
    
    pathSegment
    class(pathSegment)
    str_title
    d=["z_thres_"];
    name = pathSegment + d + str_title
    class(name)
    name = convertStringsToChars(name);
    class(name)
    saveas(fig, fullfile(filepath1, [name,'_TF_area.jpeg']));
    savefig(fig, fullfile(filepath1, [name,'_TF_area.fig']));
    
%    saveas(fig, fullfile(filepath1, [pathSegment,str_title,'_TF_area.jper']));
%    savefig(fig, fullfile(filepath1, [pathSegment,str_title,'_TF_area.fig']));
    
    area = size(query_in_area, 1)
    
    in = inpolygon(query(:, 1),query(:, 2),p.Vertices(:,1),p.Vertices(:,2)); % which coordinates of the full grid are inside the area
    tmag_xyz = sqrt(force(:,1).^2 + force(:,2).^2 + force(:,3).^2 ); %xyz
    tmag_xy = sqrt(force(:,1).^2 + force(:,2).^2 ); %xy
    [TF_mag_xyz] = griddata(V(:,1),V(:,2),tmag_xyz,query(:, 1),query(:, 2)); % interpolate the force values to the generated full grid
    [TF_mag_xy] = griddata(V(:,1),V(:,2),tmag_xy,query(:, 1),query(:, 2));
    [TFx] = griddata(V(:,1),V(:,2),force(:,1),query(:, 1),query(:, 2)); 
    [TFy] = griddata(V(:,1),V(:,2),force(:,2),query(:, 1),query(:, 2));
    [TFz] = griddata(V(:,1),V(:,2),force(:,3),query(:, 1),query(:, 2));
    
    i = 1;
%     data(i).number = i; data(i).name = list(i).name; 
    %data(i).TFxy_mean = mean(TF(:)); data(i).TFxy_max = max(TF(:)); data(i).TFxy_min  = min(TF(:));
    data(i).area = area;
    data(i).TF_in_xyz_mean = mean(TF_mag_xyz(in)); data(i).TF_in_xyz_max = max(TF_mag_xyz(in)); data(i).TF_in_xyz_min = min(TF_mag_xyz(in));
    data(i).TF_in_xy_mean = mean(TF_mag_xy(in)); data(i).TF_in_xy_max = max(TF_mag_xy(in)); data(i).TF_in_xy_min = min(TF_mag_xy(in));
    data(i).TF_in_x_mean = mean(TFx(in)); data(i).TF_in_x_max = max(TFx(in)); data(i).TF_in_x_min = min(TFx(in));
    data(i).TF_in_y_mean = mean(TFy(in)); data(i).TF_in_y_max = max(TFy(in)); data(i).TF_in_y_min = min(TFy(in));
    data(i).TF_in_z_mean = mean(TFz(in)); data(i).TF_in_z_max = max(TFz(in)); data(i).TF_in_z_min = min(TFz(in));
    writetable(struct2table(data), fullfile(filepath1, [pathSegment,'_TF_area_affected_data.xlsx']));
    
%     mean(TF); % gives NaN, not all coordinates get a force value
%     mean(TF_mag_xyz(in)) % this should be correct - gives mean traction force over the affected area
  %%  
end


function [res] = ToIdx(x, y)
    res = x*10000 + y;
end
