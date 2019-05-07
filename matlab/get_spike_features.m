function features = get_spike_features(spikes)
%% features = get_spike_features(spikes)
%
% Calculates:
% principal components
% peak valley distance
% peak height
% valley depth
% max nonlinear energy


n_spikes = length(spikes);
spike_mat = horzcat(spikes(:).GetVoltageVector())';
[eigenvectors, eigenvalues, latent] = pca(spike_mat);

features.principal_components = eigenvectors(:, 1:3);
features.pc_scores = eigenvalues(:, 1:3);

[peak_height, max_loc] = max(spike_mat, [], 2);
[valley_depth, min_loc] = min(spike_mat, [], 2);
features.peak_valley_distance = max_loc - min_loc;
features.peak_height = peak_height;
features.valley_depth = valley_depth;

max_nonlinear_energy = zeros(n_spikes, 1);
for i = 1:n_spikes
    max_nonlinear_energy(i) = max(neo(spike_mat(i, :)));
end

features.max_nonlinear_energy = max_nonlinear_energy;
