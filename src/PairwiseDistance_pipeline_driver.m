% load config info
config_info_file_path = '/home/braden/Desktop/MSKCC/PairwiseDistance_pipeline/data/config/Config.csv';
fid = fopen(config_info_file_path);
if fid < 0
    error(['could not open file: ' config_info_file_path]);
end

resource_location = '';
nuclei_file = '';
start_time = -1;
end_time = -1;

while ~feof(fid)
    % get the first line out of the way
    line = fgetl(fid);
   
    line = fgetl(fid);
    tokens = regexp(line, ',', 'split');
    
    if size(tokens, 2) == 4
        resource_location = tokens{1, 1};
        nuclei_file = tokens{1, 2};
        start_time = str2num(tokens{1, 3});
        end_time = str2num(tokens{1, 4});
    end
end

% read the nuclei file data
[embinfo, errors] = loadEmbryo_unzipped(nuclei_file, start_time, end_time);

% read the pharynx cell names
pharynx_cell_names = {};
fid = fopen(resource_location);
if fid < 0
    error(['could not open file: ' resource_location]);
end

while ~feof(fid)
    line = fgetl(fid);
    lineage_name = line(1:strfind(line, ' (')-1);
    pharynx_cell_names{end+1} = lineage_name; % consider preallocating this size and then shrinking it after loading
end
y = 0;

% create matrix with pharynx names in y, all nuc names in x
names_mat = cell(size(pharynx_cell_names, 2), size(embinfo(start_time).cellnames, 1));

% fill in names

% create a parallel matrix to hold pharynx nuc positions in y, all nuc
% positions in x
positions_mat = cell(size(pharynx_cell_names, 2), size(embinfo(start_time).cellnames, 1), 3);

% fill in the positions matrix using the nuclei data

% compute distance at every cell in the matrix

% set distance threshold

% iterate over the matrix, and add every cell whose distance is below the
% distance threshold to a list

% write the list to file
