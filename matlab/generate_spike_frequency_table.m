function generate_spike_frequency_table(mat_path, output_path, varargin)
%% generate_spike_frequency_table(mat_path, output_path, [options])
%
% Generates a table of spike frequency for each unit specified in the mat file at mat_path
%  see process_spk_files for details on generating this mat file
%
% OPTIONS
%
% bin_size - size of the time bin (in seconds) to use when counting spikes. default = 60 seconds

% anon fcn to test for files
is_file = @(fp) exist(fp, 'file');

parser = inputParser();
parser.addRequired('mat_path', is_file);
parser.addRequired('output_path');
parser.addParameter('bin_size', 60, @isnumeric);
parser.parse(mat_path, output_path, varargin{:});

mat_data = load( ...
    mat_path, ...
    'electrode_containers', ...
    'final_spike_time', ...
    'recording_start_time', ...
    'stim_times' ...
);

bin_size = parser.Results.bin_size;
electrode_containers = mat_data.electrode_containers;
final_spike_time = mat_data.final_spike_time;
recording_start_time = min(mat_data.stim_times);
%recording_start_time = mat_data.recording_start_time;
%recording_start_time = min(recording_start_time);

% only work with the containers that actually have data
containers_with_data = electrode_containers([electrode_containers(:).contains_data]);
num_units = sum([containers_with_data(:).n_clusters]);
num_bins = floor((final_spike_time - recording_start_time)/seconds(bin_size));

frequency_mat = zeros([num_bins, num_units]);

curr_unit = 1;
unit_names = {};
for curr_container = containers_with_data(:)'
    unit_names = [unit_names, curr_container.get_unit_names()];
    for iClust = 1:curr_container.n_clusters
        unit_spike_times = curr_container.spike_times( ...
            curr_container.class_no{curr_container.n_clusters} == iClust ...
        );
        frequency_mat(:, curr_unit) = generate_frequency_timecourse( ...
            unit_spike_times, ...
            'start_time', recording_start_time, ...
            'end_time', final_spike_time, ...
            'bin_size', bin_size ...
        );
        curr_unit = curr_unit + 1;
    end
end

%save('backup_mat_2.mat', 'frequency_mat', 'unit_names');
% Make sure there are no duplicate units
[unique_units, ia] = unique(unit_names);
if length(unique_units) ~= length(unit_names)
    fm2 = frequency_mat(:, ia);
    frequency_mat = fm2;
    unit_names = unique_units;
end
        
spike_table = array2table(frequency_mat, 'VariableNames', unit_names);
spike_table.time = [recording_start_time + seconds(bin_size):seconds(bin_size):final_spike_time]';
writetable(spike_table, output_path);

% Save stim times if they exist
if isfield(mat_data, 'stim_times')
    stim_times = table(mat_data.stim_times', second(mat_data.stim_times)', 'VariableNames', {'date_time', 'seconds'});
    writetable(stim_times, 'stim_times.csv');
end

%% if encounter duplicate unit name issue:
% Put breakpoint at 61: spike_table = array2table...
%unique_units = unique(unit_names);
%[unique_units, ia, ic] = unique(unit_names);
%fm2 = frequency_mat(:, ia);
%frequency_mat = fm2;
%unit_names = unique_units;
% Continue to save table.
% 
