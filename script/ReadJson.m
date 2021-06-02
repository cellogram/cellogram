function [V, D, Tv] = ReadJson(filename)

    json = fileread(filename);
    cellogram = jsondecode(json);
    clear json

    if isfield(cellogram,'analysis')
        data_fields = fieldnames(cellogram.analysis);
        for i = 1:size(data_fields)
            value = cellogram.analysis.(data_fields{i});
            if isfield(value,'data') && isfield(value,'rows')
                tmp = reshape(value.data,value.rows,value.cols);
                cellogram.analysis.(data_fields{i}) = tmp;
            end
        end
    end

    F = cellogram.analysis.F + 1;
    V = cellogram.analysis.V;
    D = cellogram.analysis.displacement;
    T = cellogram.analysis.traction_forces;
    
    % filter: only want vertices on the top surface
    Nv = size(V, 1);
    Nf = size(F, 1);
    filter = V(:, 3) > -1e-10;
    num_del = zeros(Nv, 1);
    num_del(1) = double(filter(1));
    for i = 2:Nv
        num_del(i) = num_del(i-1) + double(filter(i));
    end
    
    % force on vertex from on triangle
    force_sum = zeros(Nv, 3);
    area_sum = zeros(Nv, 1);
    for i = 1:Nf  % i-th triangle
        v1 = V(F(i, 2), :) - V(F(i, 1), :);
        v2 = V(F(i, 3), :) - V(F(i, 1), :);
        tri_area = norm(cross(v1, v2)) / 2.0;
        for j = 1:3  % j-th vertex
            if (~filter(F(i, j))), continue;end
            force_sum(F(i, j), :) = force_sum(F(i, j), :) + T(i, :) .* tri_area;
            area_sum(F(i, j)) = area_sum(F(i, j)) + tri_area;
        end
    end
    
    % output
    V = V(filter, :);
    D = D(filter, :);
    Tv = force_sum(filter, :) ./ area_sum(filter);
end