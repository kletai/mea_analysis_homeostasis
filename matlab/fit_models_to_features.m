function models = fit_models_to_features(features, varargin)
%% models = fit_models_to_features(feature, cluster_mode)
%
% Fits Gaussian Mixture Models to the specified features in the features struct
% 
% REQUIRED
% features: struct containing features to fit models to. see get_spike_features.m
%
% OPTIONAL
% cluster_mode: features to cluster on, default: 'pca' (currently no other options)
%
% PARAMS
% k_min: minimum number of clusters to compute, default = 1
% k_max: maximum number of clusters to compute, default = 5


%% parse inputs
parser = configure_parser();
parser.parse(features, varargin{:});
args = parser.Results;

switch args.cluster_mode
    case 'pca'
        fit_gm_ = @(k) fit_gm(features.pc_scores, k);
        models = arrayfun(fit_gm_, [args.k_min:args.k_max], 'UniformOutput', false);
    otherwise
        error(['Cluster Mode "', cluster_mode, '" not recognized!']);
end

if any(contains_numbers(models))
    models = models(1:find(contains_numbers(models), 1, 'first') - 1);
end


function parser = configure_parser()
    parser = inputParser;
    parser.addRequired('features');
    parser.addOptional('cluster_mode', 'pca', @ischar);
    parser.addParameter('k_max', 5);
    parser.addParameter('k_min', 1);


function model = fit_gm(X, k)
    try
        model = fitgmdist(X, k, 'Options', statset('MaxIter', 500));
    catch
        model = NaN;
    end


function result = contains_numbers(C)
    result = cellfun(@isnumeric, C);
