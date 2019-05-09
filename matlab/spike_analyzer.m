function varargout = spike_analyzer(varargin)
% SPIKE_ANALYZER MATLAB code for spike_analyzer.fig
%      SPIKE_ANALYZER, by itself, creates a new SPIKE_ANALYZER or raises the existing
%      singleton*.
%
%      H = SPIKE_ANALYZER returns the handle to a new SPIKE_ANALYZER or the handle to
%      the existing singleton*.
%
%      SPIKE_ANALYZER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SPIKE_ANALYZER.M with the given input arguments.
%
%      SPIKE_ANALYZER('Property','Value',...) creates a new SPIKE_ANALYZER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before spike_analyzer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to spike_analyzer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help spike_analyzer

% Last Modified by GUIDE v2.5 17-Jul-2017 16:34:36

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @spike_analyzer_OpeningFcn, ...
                   'gui_OutputFcn',  @spike_analyzer_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before spike_analyzer is made visible.
function spike_analyzer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to spike_analyzer (see VARARGIN)

% Choose default command line output for spike_analyzer
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes spike_analyzer wait for user response (see UIRESUME)
% uiwait(handles.figure1);
handles.plot_spike_means = true;
handles = load_data(handles);
handles = analysis_loop(handles);


% --- Outputs from this function are returned to the command line.
function varargout = spike_analyzer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

%% -------------- DATA IO -------------- %%

function handles = load_data(handles)
%% handles = load_data(handles)
%
%  Loads electrode_containers from a mat file and initializes the current container as the first
    handles.mat_path = get_mat_path();
    handles.spike_paths = get_spike_paths();
    handles.axis_loader = AxisLoader(handles.spike_paths);
    final_spike_time = load(handles.mat_path, 'final_spike_time');
    recording_start_time = load(handles.mat_path, 'recording_start_time');
    stim_times = load(handles.mat_path, 'stim_times');
    data_struct = load(handles.mat_path, 'electrode_containers');
    handles.electrode_containers = data_struct.electrode_containers;
    handles.final_spike_time = final_spike_time.final_spike_time;
    handles.recording_start_time = recording_start_time.recording_start_time;
    if isfield(stim_times, 'stim_times')
        handles.stim_times = stim_times.stim_times;
    else
        handles.stim_times = [];
    end
    handles.curr_index = 1;
    handles = load_curr_container(handles);
    if ~handles.curr_container.contains_data
        disp_skip_msg(handles);
        handles = next_electrode(handles);
    end

function mat_path = get_mat_path()
    [mat_filename, mat_dirname] = uigetfile('*.mat', 'Select matfile to load')
    mat_path = fullfile(mat_dirname, mat_filename);

function spike_paths = get_spike_paths()
    [spike_filenames, spike_dirname] = uigetfile('*.spk', 'Select spkfile(s) to load', 'MultiSelect', 'on');
    % uigetfile returns a cell if multiple files are selected but a string if only one file is selected
    if iscell(spike_filenames) 
        spike_filenames = spike_filenames;
    else
        spike_filenames = {spike_filenames};
    end
    full_spike_path = @(spike_filename) fullfile(spike_dirname, spike_filename);
    spike_paths = cellfun(full_spike_path, spike_filenames, 'UniformOutput', false)

function handles = load_curr_container(handles)
    electrode_container = handles.electrode_containers(handles.curr_index);
    handles.curr_container = electrode_container.create_spike_container(handles.axis_loader);
    n_spikes = handles.curr_container.get_number_of_spikes();
    n_sample_spikes = 1000;
    if n_spikes > n_sample_spikes 
        handles.sample_spikes = datasample(1:n_spikes, n_sample_spikes, 'Replace', false);
    else
        handles.sample_spikes = 1:n_spikes;
    end

function handles = update_curr_container(handles)
    handles.electrode_containers(handles.curr_index) = handles.curr_container.get_electrode_container();

function save_data(handles)
%% save_data(handles)
%
%
    electrode_containers = handles.electrode_containers;
    final_spike_time = handles.final_spike_time
    recording_start_time = handles.recording_start_time
    stim_times = handles.stim_times;
    save(handles.mat_path, 'electrode_containers', 'final_spike_time', 'recording_start_time', 'stim_times', '-v7.3');

%% --------------- DISPLAY ----------------%%

function handles = refresh_display(handles)
%% refresh_display(handles)
%
% 
    axes(handles.main_axes)
    handles = plot_features(handles);
    axes(handles.secondary_axes)
    if handles.plot_spike_means
        handles = plot_spike_means(handles);
    else
        handles = plot_spikes(handles);
    end
    handles = refresh_colors(handles);

function handles = plot_features(handles)
%% plot_features(container)
%
% Creates a 3d scatter plot of the data in pc space
    features = handles.curr_container.features.pc_scores(handles.sample_spikes, :);
    handles.scatter_handle = scatter3(features(:, 1), features(:, 2), features(:, 3), 'filled');
    labels = arrayfun(@num2str, 1:handles.curr_container.n_clusters, 'UniformOutput', false);
    legend(labels)

function handles = plot_spikes(handles)
%% plot_spikes(container)
%
%
    handles.spike_handles = cell(size(handles.sample_spikes));
    spikes = handles.curr_container.get_spikes();
    spike_no_all = [];
    for i = 1:length(handles.sample_spikes)
        spike_no = handles.sample_spikes(i);
        handles.spike_handles{i} = plot(spikes(spike_no, :));
        spike_no_all = [spike_no_all, spike_no];
        hold on
    end
    hold off

function handles = plot_spike_means(handles)
    handles.spike_mean_handles = cell(size(handles.curr_container.n_clusters, 1));
    mean_spikes = handles.curr_container.mean_waveforms{handles.curr_container.n_clusters};
    for iCluster = 1:handles.curr_container.n_clusters
        handles.spike_mean_handles{iCluster} = plot(mean_spikes(iCluster, :), 'LineWidth', 2);
        hold on
    end
    hold off

function handles = refresh_colors(handles)
%% handles = refresh_colors(handles)
%
%
    n_samples = size(handles.sample_spikes, 1);
    colors = lines(handles.curr_container.n_clusters);
    c_data = zeros(n_samples, 3);
    sample_classes = handles.curr_container.class_no{handles.curr_container.n_clusters}(handles.sample_spikes);
    no_class_color = [0.7, 0.7, 0.7];
    for i = 1:handles.curr_container.n_clusters
        if ~handles.curr_container.valid_class(i)
            colors(i, :) = no_class_color;
        end
    end

    if handles.curr_container.n_clusters == 1
        nspikes = size(handles.scatter_handle.XData,2);
        blue = [0, 0, 1];
        light = [203, 192, 255]/255;
        c_gra = [linspace(blue(1),light(1),nspikes)', linspace(blue(2),light(2),nspikes)', linspace(blue(3),light(3),nspikes)'];
        [~, sorted_indices] = sort(handles.curr_container.spike_times(handles.sample_spikes));
        c_data = c_gra(sorted_indices,:);
    else
        for i = 1:length(sample_classes)
            c_data(i, :) = colors(sample_classes(i), :);
        end
    end

    if handles.plot_spike_means
        for i = 1:handles.curr_container.n_clusters
            handles.spike_mean_handles{i}.Color = colors(i, :);
        end
    else
        for i = 1:length(sample_classes)
            handles.spike_handles{i}.Color = colors(sample_classes(i), :);
        end
    end

    handles.scatter_handle.CData = c_data;

%% ----------------- MISC ---------------- %%
function disp_progress_msg(handles)
    disp([ ...
       'Analyzing electrode ', ...
       num2str(handles.curr_index), ...
       ' of ', ...
       num2str(numel(handles.electrode_containers))...
    ]);

function disp_skip_msg(handles)
    disp([
        'skipping electrode ', ...
        get_well_string( ... 
            handles.electrode_containers(handles.curr_index).spike_index(1), ...
            handles.electrode_containers(handles.curr_index).spike_index(2) ...
        ), ...
        '-', ...
        get_electrode_string( ...
            handles.electrode_containers(handles.curr_index).spike_index(3), ...
            handles.electrode_containers(handles.curr_index).spike_index(4) ...
        ) ...
    ]);


%% -------------- MAIN LOOP -------------- %%

function command = get_command(handles)
%% command = get_command(handles)
%
%  Retrieves the next command from the input prompt
    axes(handles.main_axes);
    [~, ~, command] = ginput(1);

function handles = set_n_clust(handles, n_clust)
%% handles = set_n_clust(handles)
%
%  Queries the user for the number of clusters to use
    max_clusters = length(handles.curr_container.cluster_model);
    if n_clust > 0 && n_clust <= max_clusters
        handles.curr_container.n_clusters = n_clust;
        handles.curr_container.valid_class = true(n_clust, 1);
    else
        disp(['Number of clusters must be between 0 and ', num2str(max_clusters), '!'])
    end

function handles = merge_clusters(handles, merging_clusters)
%
%   Merges clusters specified by user - useful when the automatic
%   clustering appears inaccurate

    clust_nums = find(merging_clusters == 1); % Cluster numbers that are merging
    if ~isempty(clust_nums)
        num_merging = length(clust_nums); % How many clusters are merging?
        new_n_clusters = handles.curr_container.n_clusters - num_merging + 1; % How many clusters will there be now?
        merge_to = clust_nums(1);
        merge_from = handles.curr_container.n_clusters; % Number of clusters before merging
        all_nums = handles.curr_container.class_no{merge_from};
        % Change the class_no of merging clusters
        for i = 2:num_merging 
            change_inds = find(all_nums == clust_nums(i));
            for j = change_inds
                all_nums(j) = merge_to;
            end
        end
        % Check that class_no doesn't skip cluster numbers now (cluster numbers
        % should go 1:new_n_clusters)
        if max(all_nums) > new_n_clusters
            kept_nos = unique(all_nums); % Finds which cluster numbers are still stored
            removed_nos = clust_nums(2:end); % Cluster numbers that were removed in merge
            over_kept = kept_nos(kept_nos > new_n_clusters); % Finds kept class_nos higher than max class_no
            under_rem = removed_nos(removed_nos <= new_n_clusters); % Finds removed class_nos that need to be replaced
            replace = find(all_nums > new_n_clusters);
            for i = 1:length(replace)
                replace_ind = replace(i);
                replace_with = find(over_kept == all_nums(replace_ind));
                all_nums(replace_ind) = under_rem(replace_with);
            end
        end    
        % Update handles, saving the old unmerged class_no in cell 6 of
        % handles.curr_container.class_no
        handles.curr_container.class_no{6} = handles.curr_container.class_no{new_n_clusters};
        handles.curr_container.class_no{new_n_clusters} = all_nums;
        handles.curr_container.n_clusters = new_n_clusters;
    end

function handles = undo_merge(handles)
%
%   Returns to original clustering determined during preprocessing
    if max(handles.curr_container.class_no{6}) == handles.curr_container.n_clusters
        handles.curr_container.class_no{handles.curr_container.n_clusters} = handles.curr_container.class_no{6};
        handles.curr_container.class_no{6} = [];
    else
        disp('No clusters were merged on this plot');
    end

function print_details(handles)
    disp(handles.curr_container)

function handles = flag_clusters(handles)
    handles.curr_container.valid_class = flag_popup(handles);

function rotate_display(handles)
    axes(handles.main_axes)
    rotate3d()
    pause
    rotate3d()

function handles = next_electrode(handles)
    if handles.curr_index < numel(handles.electrode_containers)
        handles.curr_index = handles.curr_index + 1;
    else
        disp('Wrapping back to first electrode!')
        handles.curr_index = 1;
    end
    if handles.electrode_containers(handles.curr_index).contains_data
        handles = load_curr_container(handles);
        disp_progress_msg(handles);
    else
        % this will overflow if all of the electrodes are empty
        disp_skip_msg(handles);
        handles = next_electrode(handles);
    end

function handles = prev_electrode(handles)
    if handles.curr_index > 1
        handles.curr_index = handles.curr_index - 1;
    else
        disp('Wrapping to last electrode!')
        handles.curr_index = numel(handles.electrode_containers);
    end
    if handles.electrode_containers(handles.curr_index).contains_data
        handles = load_curr_container(handles);
        disp_progress_msg(handles);
    else
        % this will overflow if all of the electrodes are empty
        disp_skip_msg(handles);
        handles = prev_electrode(handles);
    end
    
function handles = choose_electrode(handles)
    chosen_index = sel_index_popup();
    if ~isempty(chosen_index) && chosen_index < numel(handles.electrode_containers) && chosen_index > 0
        handles.curr_index = chosen_index;
        if handles.electrode_containers(handles.curr_index).contains_data
            handles = load_curr_container(handles);
            disp_progress_msg(handles);
        else
            % this will overflow if all of the electrodes are empty
            disp_skip_msg(handles);
            handles = next_electrode(handles);
        end
    else
        disp('Invalid electrode.')
    end

function handles = analysis_loop(handles)
    keep_looping = true;
    handles = refresh_display(handles);
    while keep_looping
        command = get_command(handles);
        switch command
            case {'1', '2', '3', '4', '5'}
                n_clust = uint8(command) - 48
                handles = set_n_clust(handles, n_clust);
                if handles.plot_spike_means
                    axes(handles.secondary_axes)
                    handles = plot_spike_means(handles);
                end
                handles = refresh_colors(handles);
            case 'm'
                % Select clusters to be merged
                merging_clusters = sel_clust_popup(handles);
                % Merge clusters
                handles = merge_clusters(handles, merging_clusters);
                handles = refresh_colors(handles);
            case 'u'
                % Undo merge
                handles = undo_merge(handles);
                handles = refresh_colors(handles);
            case 'd'
                print_details(handles);
            case 'f'
                handles = flag_clusters(handles);
                handles = refresh_colors(handles);
            case 'r'
                rotate_display(handles);
            case 'v'
                handles.plot_spike_means = true;
                axes(handles.secondary_axes);
                handles = plot_spike_means(handles);
                handles = refresh_colors(handles);
            case 't'
                handles.plot_spike_means = false;
                axes(handles.secondary_axes);
                handles = plot_spikes(handles);
                handles = refresh_colors(handles);
            case 'x'
                handles = update_curr_container(handles);
                handles = next_electrode(handles);
                handles = refresh_display(handles);
            case 'z'
                handles = update_curr_container(handles);
                handles = prev_electrode(handles);
                handles = refresh_display(handles);
            case 'g'
                handles = update_curr_container(handles);
                handles = choose_electrode(handles);
                handles = refresh_display(handles);                
            case 's'
                handles = update_curr_container(handles);
                save_data(handles)
            otherwise
                disp(['Command: "', command, '" not recognized!'])
        end
    end