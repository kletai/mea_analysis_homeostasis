function pca_spikes(data_path, electrode_index)

axis_file = AxisFile(data_path);
data = axis_file.DataSets.LoadData();

electrode_waveforms = data(electrode_index(1), electrode_index(2), electrode_index(3), electrode_index(4));
waveform_mat = single(horzcat(electrode_waveforms{:}.Data)');

figure
for i = 1:size(waveform_mat, 1)
    plot(waveform_mat(i, :))
    hold on
end

[coeff, score, latent] = pca(waveform_mat);

figure
plot(latent)
title('explained variance')

figure
scatter3(score(:, 1), score(:, 2), score(:, 3))
xlabel('PC1')
ylabel('PC2')

center = mean(score(:, [1,2,3]), 1);
hold on
scatter3(center(1), center(2), center(3), '*r')
hold off

figure
for i = 1:3
    plot(coeff(i, :))
    hold on
end
legend({'1', '2', '3'})

figure
plot(center(1)*coeff(1, :) + center(2)*coeff(2, :) + center(3)*coeff(3, :))


% Fit a GMM
k = 5
%gmm = fitgmdist(score(:, 1:3), 3)
[cluster_inds, cluster_mus] = kmeans(score(:, 1:3), k);
figure
for i = 1:k
    matching_points = cluster_inds == i;
    scatter3(score(matching_points, 1), score(matching_points, 2), score(matching_points, 3))
    hold on
end

figure
for i = 1:k
    plot(cluster_mus(i, 1)*coeff(1, :) + cluster_mus(i, 2)*coeff(2, :) + cluster_mus(i, 3)*coeff(3, :))
    hold on
end
