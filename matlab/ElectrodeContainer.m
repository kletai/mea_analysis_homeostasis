classdef ElectrodeContainer < handle
    properties
        well
        electrode
        spike_index
        spike_times
        class_no
        valid_class
        features
        cluster_model
        n_clusters
        contains_data
        mean_waveforms
        % If you add a property here, you must add it to the constructor calls in ElectrodeContainer.create_spike_container() and SpikeContainer.get_electrode_container()
    end
    
    methods
        function obj = ElectrodeContainer(varargin)
            %% ElectrodeContainer([spike_index], varargin)
            %
            % Container for data associated with an electrode
            %   
            % PARAMETERS
            %
            % spike_index - contains the [well_row, well_col, electrode_col, electrode_row] index 
            %   corresponding to the electrode that the  ElectrodeContainer stores data for (TODO - rename this to electrode_index)
            %
            % spike_times - contains the timing of each spike (as a matlab datetime object) that occured on the electrode
            %
            % contains_data - true for electrodes containing a sufficient number of spikes for preprocessing. 
            %   Electrodes with fewer than 20 spikes are ignored during preprocessing. 
            %
            % valid - n_clusters x 1 boolean array set using the 'flag' command on spike_analyzer. 
            %   false for clusters the user has determined to be unusable
            %
            % class_no - num_clustering_results x 1 cell array containing the n_spikes x 1 classification arrays 
            %   for each clustering. For example, to access the class of the 5th spike for a 3 cluster model, 
            %   you can do this: spike_class = ElectrodeContainer.class_no{3}(5). 
            %   The classification in this case would be a number from 1 - 3.
            %
            % features - struct containing spike features. features.pc_scores contains eigenvalues for 
            %   the first 3 principal component.s this is the only feature that I use right now, 
            %   although there are others computed during preprocessing. See get_spike_features.m
            %
            % cluster_model - num_clustering_results x 1 cell array containing model objects for each clustering. 
            %   currently contains matlab gmdistribution objects
            %
            % n_clusters - the number of clusters determined to be present on an electrode. 
            %   Use this as an index to the cluster_model and class_no fields to retrieve the 
            %   appropriate cluster model object and spike classifications respectively.
            %
            % mean_waveforms - num_clustering_results x 1 cell array containing num_clusters x spike_dimensionality 
            %   waveform matricies containing the mean waveforms for each cluster for each clustering result.

            parser = inputParser();
            parser.addOptional('spike_index', [], @isnumeric);
            parser.addOptional('spike_times', []);
            parser.addParameter('contains_data', true);
            parser.addParameter('valid', true);
            parser.addParameter('class_no', 0);
            parser.addParameter('features', false);
            parser.addParameter('cluster_model', false);
            parser.addParameter('n_clusters', 1);
            parser.addParameter('mean_waveforms', []);
            
            parser.parse(varargin{:});
            obj.spike_index = parser.Results.spike_index;
            obj.spike_times = parser.Results.spike_times;
            obj.contains_data = parser.Results.contains_data;
            obj.class_no = parser.Results.class_no;
            obj.valid_class = parser.Results.valid;
            obj.features = parser.Results.features;
            obj.cluster_model = parser.Results.cluster_model;
            obj.n_clusters = parser.Results.n_clusters;
            obj.mean_waveforms = parser.Results.mean_waveforms;
        end

        function names = get_unit_names(obj)
            %% names = get_unit_names(obj)
            %
            % returns a cell array of unit names of the form [well_string, electrode_string, unit_number]
            %
            % ex. Unit 1 on electrode A1-12 would be 'A1121'

            well_ele_name = [ ...
                get_well_string(obj.spike_index(1), obj.spike_index(2)), ...
                get_electrode_string(obj.spike_index(3), obj.spike_index(4)) ...
            ];
            name_generator = @(n) sprintf([well_ele_name, '%d'], n);
            names = arrayfun(name_generator, 1:obj.n_clusters, 'UniformOutput', false);
        end

        function spk_container = create_spike_container(obj, spikes)
            spk_container = SpikeContainer( ...
                spikes, ...
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

        function n_spikes = get_number_of_spikes(obj)
            if obj.contains_data
                n_spikes = size(obj.features.pc_scores, 1);
            else
                n_spikes = 0;
            end
        end
    end
end
