classdef SpikeContainer < ElectrodeContainer
    %% SpikeContainer < Electrode Container
    % This handles loading of the spike data corresponding to an ElectrodeContainer
    properties (Access = public)
        data_loaded
    end
    properties (Access = private)
        axis_loader
        spikes
    end
    methods
        function obj = SpikeContainer(axis_loader, varargin)
            %%  obj = SpikeContainer(spikes, varargin)
            %
            % spikes - n x d containing spike waveforms, where n is the number of spikes, 
            %  and d is the dimensionality of each spike

            %% Call superclass constructor
            obj = obj@ElectrodeContainer(varargin{:})

            %%
            parser = inputParser();
            parser.addRequired('axis_loader', @(al) isa(al, 'AxisLoader'))
            parser.parse(axis_loader);
            obj.axis_loader = parser.Results.axis_loader;
            obj.data_loaded = false;
            obj.spikes = [];
        end

        function electrode_container = get_electrode_container(obj)
            electrode_container = ElectrodeContainer( ...
                obj.spike_index, ...
                obj.spike_times, ...
                'contains_data', obj.contains_data, ...
                'valid', obj.valid_class, ...
                'class_no', obj.class_no, ...
                'features', obj.features, ...
                'cluster_model', obj.cluster_model, ...
                'n_clusters', obj.n_clusters, ...
                'mean_waveforms', obj.mean_waveforms ...
            );
        end

        function load_spike_data(obj)
            spk_data = obj.axis_loader.load_data_set_from_index(obj.spike_index);
            obj.spikes = horzcat(spk_data(:).GetVoltageVector())';
            obj.data_loaded = true;
        end

        function spikes = get_spikes(obj)
            if ~obj.data_loaded
                obj.load_spike_data();
            end
            spikes = obj.spikes;
        end

        function avg_waveforms = get_average_waveforms(obj)
            if ~obj.data_loaded
                obj.load_spike_data();
            end
            avg_waveforms = zeros([obj.n_clusters, size(obj.spikes, 2)])
            for i = 1:obj.n_clusters
                avg_waveforms(i, :) = mean(obj.spikes(obj.class_no == i, :), 1);
            end
        end
    end
end
