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
names_mat = cell(size(pharynx_cell_names, 2)+1, size(embinfo(start_time).cellnames, 1)+1);

% fill in names

% create a parallel matrix to hold pharynx nuc positions in y, all nuc
% positions in x
positions_mat = cell(size(pharynx_cell_names, 2)+1, size(embinfo(start_time).cellnames, 1)+1, 3);

% fill in the names and positions in the matrices using the nuclei data
% access the nuclei info at the time point
cell_data_atT = embinfo(start_time).celldata;
cell_names = embinfo(start_time).cellnames;

% make sure the arrays are parallel in the number of rows
mat_x_itr = 2;
mat_y_itr = 2;
if size(cell_data_atT, 1) == size(cell_names, 1)
    % iterate over all the cell names
    for i=1:size(cell_names, 1)
        name = cell_names{i, 1};
        cell_data = cell_data_atT(i, :);
        x = cell_data(4);
        y = cell_data(5);
        z = cell_data(6);
        diam = cell_data(7);
        
        % add the name to the x dimension of names_mat and the positional
        % info to the x dimension of the positions mat
        names_mat{1, mat_x_itr} = name;
        positions_mat(1, mat_x_itr, :) = {x, y, z};
        mat_x_itr = mat_x_itr + 1;
           
        % check if the lineage name is part of the pharynx by checking the
        % pharynx names file
        idx = find(strcmp(pharynx_cell_names, name), 1);
        if ~isempty(idx)
            % add the name to the y dimension of the names_mat and the
            % positional info to the y dimension of the positions mat
            names_mat{mat_y_itr, 1} = name;
            positions_mat(mat_y_itr, 1, :) = {x, y, z};
            mat_y_itr = mat_y_itr + 1;
        end
    end
end

% set distance threshold
dist_threshold = 20;

% set up matrix to record distances under threshold:
% non-pharynx nuc - pharynx nuc - distance
pharynx_neighbors = cell(50, 3);
pharynx_neighbors_itr = 1;

% compute distance at every cell in the matrix
% record every distance calculation in names_mat table, and save those
% under the threshold in 
% iterate over rows
for y=2:size(names_mat, 1)
   % iterate over columns
   for x=2:size(names_mat, 2)
       
      % compute the distance between the two nuclei 
      x1 = positions_mat(y, 1, 1);
      y1 = positions_mat(y, 1, 2);
      z1 = positions_mat(y, 1, 3);
      
      x2 = positions_mat(1, x, 1);
      y2 = positions_mat(1, x, 2);
      z2 = positions_mat(1, x, 3);
       
      pointsXdimensions = [x1 y1 z1; x2 y2 z2];
      d = pdist(cell2mat(pointsXdimensions));
      
      % insert distance in names_mat
      names_mat{y, x} = d;
      
      % check if under distance threshold
      if d < dist_threshold
      
        % make sure not distance between two pharynx nuclei
        idx = find(strcmp(pharynx_cell_names, names_mat{1, x}), 1);
        if isempty(idx)
            
            % make sure not already recorded nuc
            idx = find(strcmp(pharynx_neighbors(:, 1), names_mat{1, x}), 1);
            if isempty(idx)
               
                % add the non-pharynx name, pharynx cell neighbor, and
                % distance to pharynx_neighbors list
                pharynx_neighbors(pharynx_neighbors_itr, :) = {names_mat{1, x}, names_mat{y, 1}, d};
                pharynx_neighbors_itr = pharynx_neighbors_itr + 1;
            end
        end
      end
   end    
end

% write the list to file
neighbors_list_fileID = fopen('/home/braden/Desktop/MSKCC/PairwiseDistance_pipeline/data/output/pharynx_neighbors_list.txt', 'w');
fprintf(neighbors_list_fileID, '%s, %s, %s','neighbor', 'phar nuc' , 'distance');
fprintf(neighbors_list_fileID, '\n');
for i=1:size(pharynx_neighbors, 1)
    
    fprintf(neighbors_list_fileID, '%s, %s, %s',pharynx_neighbors{i, 1}, pharynx_neighbors{i, 2}, num2str(pharynx_neighbors{i, 3}));
    fprintf(neighbors_list_fileID, '\n');
end
fclose(neighbors_list_fileID);

% write the names_mat to file


% format the list of neighbors into a wormguides readable URL
url_start = 'http://scene.wormguides.org/wormguides/testurlscript?/set/e2-n$+#ff4d1a4d/e1-n$+#ff4d1a4d/e3-n$+#ff4d1a4d/mc1-n$+#ff33cc53/mc2-n$+#ff33cc53/mc3-n$+#ff33cc53/g1a-n$+#ff33cc53/g1p-n$+#ff33cc53/g2r-n$+#ff33cc53/g2l-n$+#ff33cc53/mcl-n$+#ff33cc53/i1-n$+#ff33cc53/i2-n$+#ff33cc53/nsm-n$+#ff33cc53/i3-n$+#ff33cc53/i4-n$+#ff33cc53/i5-n$+#ff33cc53/i6-n$+#ff33cc53/m3-n$+#ff33cc53/m2-n$+#ff33cc53/m1-n$+#ff33cc53/m4-n$+#ff33cc53/m5-n$+#ff33cc53/m6-n$+#ff33cc53/m7-n$+#ff33cc53/vpi-n$+#ffffffb3/buccal-d$+#ffffffb3/nerve_ring_right-M+#ffffffff/nerve_ring_left-M+#ffffffff/int-n$+#ffffffff/E-s<$+#ffffff4d/MSaa-s<$+#ffffff4d/MSpaa-s<$+#ffffff4d/MSpapa-s<$+#ffffff4d/ABalpaa-s<$+#ffffff4d/ABalpapp-s<$+#ffffff4d/ABaraa-s<$+#ffffff4d/ABarapa-s<$+#ffffff4d/Pharynx-M+#8199ff00/'
url_end = 'view/time=251/rX=-44.0/rY=14.0/rZ=0.0/tX=2.75/tY=-1.75/scale=2.0/dim=0.25/browser/'
nuc_rule_format = '-s$+#ffe64d4d/';

url = url_start;
for i=1:size(pharynx_neighbors, 1)
    % access the lineage name
    nuc_name = pharynx_neighbors{i, 1};
    
    % make sure not an empty cell
    if isempty(nuc_name)
        continue
    end
    
    % format the lineage name with the rule
    new_rule = strcat(nuc_name, nuc_rule_format);
    
    % concatenate the rule with the rest of the url
    url = strcat(url, new_rule);
end

% append the url with the end string
url = strcat(url, url_end);

% write the URL to a file
url_fileID = fopen('/home/braden/Desktop/MSKCC/PairwiseDistance_pipeline/data/output/pharynx_neighbors_wormguides_visualization_url.txt', 'w');
fprintf(url_fileID, url);
fclose(url_fileID);



