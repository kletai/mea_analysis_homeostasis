classdef AxisLoader < handle
    properties
        data_paths
        file_objs
        channels
        channel_ind
        num_channels
        recording_start_time
    end
    
    methods
        function obj = AxisLoader(data_paths, varargin)
            %% obj = AxisLoader(data_paths)
            %
            % An object for loading data from multiple
            parser = inputParser();
            parser.addRequired('data_paths', @iscell);
            parser.parse(data_paths, varargin{:});
            results = parser.Results;
            
            obj.data_paths = results.data_paths;
            obj.file_objs = cellfun(@AxisFile, results.data_paths, 'UniformOutput', false);
            
            obj.assert_consistent_plate_type();
            obj.channels = obj.file_objs{1}.DataSets.ChannelArray.Channels;
            obj.channel_ind = 0;
            obj.num_channels = numel(obj.channels);
            obj.recording_start_time = [];%datetime('3000-10-06 00:00:00'); % hbd
            for i = 1:numel(obj.file_objs)
                file_start_time = spike_time_to_datetime(obj.file_objs{i}.DataSets.Header.FileStartTime);
                obj.recording_start_time = [obj.recording_start_time, file_start_time];
            end
        end
        
        function assert_consistent_plate_type(obj)
            %% AxisLoader.assert_consistent_plate_type()
            %
            % raises an error if the plate types for each AxisFile in file_objs are not the same
            for i = 1:length(obj.file_objs)
                if obj.file_objs{i}.DataSets.ChannelArray.PlateType ~= obj.file_objs{1}.DataSets.ChannelArray.PlateType
                    error('spk files must be recorded from the same plate type!')
                end
            end
        end
        
        function [data_set, next_channel] = load_next_data_set(obj)
            %%  [data_set, loaded_channel] = AxisLoader.load_next_dataset(obj)
            %
            % loads the next data set from the axis files
            %  note that data sets are not necesarily loaded in a
            if obj.channel_ind < obj.num_channels
                obj.channel_ind = obj.channel_ind + 1;
            else
                error('All channels have been loaded!')
            end
            next_channel = obj.channels(obj.channel_ind);
            data_set = obj.load_data_set(next_channel);
        end
        
        function data_set = load_data_set_from_index(obj, index)
            %% data_set = load_data_set_from_index(obj, index)
            %
            % loads data from a vectorized index, rather than a channel object
            channel = struct( ...
                'WellRow', index(1), ...
                'WellColumn', index(2), ...
                'ElectrodeColumn', index(3), ...
                'ElectrodeRow', index(4) ...
            );
            data_set = obj.load_data_set(channel);
        end
        
        function data_set = load_data_set(obj, channel)
            %% data_set = AxisLoader.load_data_set(obj, channel)
            %
            % loads the data_set at the channel specified by channel
            data_cell = obj.load_data_cell(channel);
            data_set = data_cell{ ...
                channel.WellRow, ...
                channel.WellColumn, ...
                channel.ElectrodeColumn, ...
                channel.ElectrodeRow ...
            };
        end
        
        function data_cell = load_data_cell(obj, channel)
            %% data_cell = AxisLoader.load_data_cell(obj, channel)
            %
            % loads an well_rows x well_cols x elec_cols x elec_rows
            %  containing the pooled data from each AxisFile in obj.file_objs
            %  for the specified channel
            well_string = get_well_string( ...
                channel.WellRow, ...
                channel.WellColumn ...
            );
            electrode_string = get_electrode_string( ...
                channel.ElectrodeColumn, ...
                channel.ElectrodeRow ...
            );
            % Even though we are only loading data from individual electrodes,
            %  a full sized cell array is returned
            transfer_cells = cell(size(obj.file_objs));
            for i = 1:numel(transfer_cells)
                transfer_cells{i} = obj.file_objs{i}.DataSets.LoadData(well_string, electrode_string);
            end
            % This concatenates the contents of the cells loaded from each AxisFile,
            %  on an electrode by electrode basis
            data_cell = cellfun(@horzcat, transfer_cells{:}, 'UniformOutput', false);
        end
        
        function cell_shape = get_cell_shape(obj)
            channels = obj.file_objs{1}.DataSets.ChannelArray.Channels;
            cell_shape = [ ...
                max([channels(:).WellRow]), ...
                max([channels(:).WellColumn]), ...
                max([channels(:).ElectrodeColumn]), ...
                max([channels(:).ElectrodeRow]) ...
            ];
        end
    end
end
