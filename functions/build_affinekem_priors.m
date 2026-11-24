function [K, composite_reference, group_for_frame, group_duration] = build_affinekem_priors( ...
    project, backproject, counts, normalization, background, durations, image_size, prior_iterations)
%BUILD_AFFINEKEM_PRIORS Construct the sparse kernel and composite references.
fprintf('Building composite-frame kernel features and affine references.\n');
groups = {1:16, 17:20, 21:24};
features = zeros(prod(image_size), numel(groups));
composite_reference = zeros(size(features));
group_for_frame = zeros(1,24);
group_duration = zeros(1,numel(groups));
feature_filter = fspecial('gaussian', [3 3], 3/(4*sqrt(2*log(2))));
for group = 1:numel(groups)
    members = groups{group};
    group_for_frame(members) = group;
    group_duration(group) = sum(durations(members));
    y = sum(counts(:,members), 2);
    r = sum(background(:,members), 2);
    n = repmat(normalization, 1, numel(members));

    kernel_prior = composite_em(project, backproject, y, sum(n, 2), r, ...
        prod(image_size), prior_iterations);
    smooth_prior = convn(reshape(kernel_prior, image_size), feature_filter, 'same');
    features(:,group) = smooth_prior(:) / group_duration(group);
    composite_reference(:,group) = composite_em(project, backproject, y, mean(n, 2), r, ...
        prod(image_size), prior_iterations);
end
feature_std = std(features, 1, 1);
assert(all(feature_std > 0), ...
    'A composite kernel feature is constant. Increase TotalCounts or change Seed.');
features = features ./ feature_std;
[neighbors, weights] = buildKernel(image_size, 'knn', 48, features, 'radial', 1);
K = buildSparseK(neighbors, weights);
end

function x = composite_em(project, backproject, counts, normalization, background, n_pixels, iterations)
% ML-EM step used to construct the composite image inputs.
x = ones(n_pixels, 1);
sensitivity = backproject(normalization);
sensitivity(sensitivity <= 0) = 1e-12;
floor_prediction = max(mean(counts) * 1e-9, 1e-12);
for iteration = 1:iterations
    predicted = max(normalization .* project(x) + background, floor_prediction);
    correction = backproject(normalization .* (counts ./ predicted));
    correction(~isfinite(correction) | correction < 0) = 0;
    x = x ./ sensitivity .* correction;
    x(~isfinite(x) | x < 0) = 0;
end
end
