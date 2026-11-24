function [counts, normalization, background, truth, durations, image_size, project, backproject] = ...
    simulate_zubal_pet(data_file, total_counts)
%SIMULATE_ZUBAL_PET Generate projections and ground truth for all 24 frames.
assert(isfile(data_file), 'Zubal data not found. Set DataFile to zubal_data.mat.');
data = load(data_file, 'model', 'TAC', 'scant', 'u');
image_size = [111 111];
assert(all(isfield(data, {'model','TAC','scant','u'})), ...
    'The Zubal MAT file must contain model, TAC, scant, and u.');
assert(isequal(size(data.model), image_size) && isequal(size(data.u), image_size) ...
    && isequal(size(data.TAC), [24 6]) && isequal(size(data.scant), [24 2]), ...
    'Expected a 111-by-111 phantom and attenuation map, with 24 frames and six TACs.');

durations = (data.scant(:,2) - data.scant(:,1)) / 60;
assert(all(isfinite(durations) & durations > 0), 'Frame durations must be positive and finite.');

fprintf('Building the 2-D PET projector and simulating 24 frames.\n');
ig = image_geom('nx', image_size(1), 'ny', image_size(2), 'fov', 33.3);
ig.mask = ig.circ(16, 16) > 0;
sg = sino_geom('par', 'nb', 249, 'na', 210, 'dr', 70/249, 'strip_width', 'dr');
G = Gtomo2_strip(sg, ig, 'single', 1);
mask = ig.mask(:);
project = @(x) double(G * x(mask));
backproject = @(y) reshape(double(ig.embed(double(G' * y))), [], 1);

labels = [0 32 64 128 196 255];
truth = zeros(prod(image_size), 24);
for tissue = 1:numel(labels)
    region = data.model(:) == labels(tissue);
    truth(region,:) = repmat((data.TAC(:,tissue) .* durations)', nnz(region), 1);
end

projection = zeros(sg.nb * sg.na, 24);
for frame = 1:24
    projection(:,frame) = project(truth(:,frame));
end
attenuation = exp(-project(data.u(:)));
signal = attenuation .* projection;
background = repmat(0.2 * mean(signal, 1), size(signal,1), 1);
expected = signal + background;
scale = total_counts / sum(expected(:));
counts = poissrnd(scale * expected);
background = scale * background;
normalization = scale * attenuation;
end
