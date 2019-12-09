clear
[file,filepath]=uigetfile({'*.json*','JSON Files'},...
    'Select Data File');

json = fileread(fullfile(filepath,file));
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

if isfield(cellogram,'hull')
    data_fields = fieldnames(cellogram.hull);
    for i = 1:size(data_fields)
        value = cellogram.hull.(data_fields{i});
        if isfield(value,'data') && isfield(value,'rows')
            tmp = reshape(value.data,value.rows,value.cols);
            cellogram.hull.(data_fields{i}) = tmp;
        end
    end
end

if isfield(cellogram,'mesh')
    data_fields = fieldnames(cellogram.mesh);
    for i = 1:size(data_fields)
        value = cellogram.mesh.(data_fields{i});
        if isfield(value,'data') && isfield(value,'rows')
            tmp = reshape(value.data,value.rows,value.cols);
            cellogram.mesh.(data_fields{i}) = tmp;
        end
    end
end

filename = strsplit(file,'.');

clear data_fields tmp value i

save(fullfile(filepath,filename{1}),'cellogram');

figure
% trimesh(cellogram.analysis.F+1,cellogram.analysis.V(:,1),cellogram.analysis.V(:,2),cellogram.analysis.V(:,3))
trimesh(cellogram.mesh.triangles+1,cellogram.mesh.points(:,1),cellogram.mesh.points(:,2),cellogram.mesh.points(:,3))
%tetramesh(cellogram.analysis.T(1:200,:)+1,cellogram.analysis.V)