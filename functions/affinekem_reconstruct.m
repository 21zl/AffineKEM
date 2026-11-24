function history = affinekem_reconstruct(project, backproject, K, counts, normalization, ...
    background, reference, lambda, record_iterations)
%AFFINEKEM_RECONSTRUCT Reconstruct one frame using a fixed composite reference.
% Inputs are column vectors. record_iterations is increasing and positive.
% project and backproject map between full image and projection vectors.

pixel_sensitivity = backproject(normalization);
pixel_sensitivity(pixel_sensitivity < 0) = 0;
support = pixel_sensitivity > 0;
kernel_sensitivity = K' * pixel_sensitivity;
kernel_sensitivity(kernel_sensitivity <= 0) = 1e-12;

alpha = double(support);
bias = zeros(size(alpha));
kernel_image = K * alpha;
kernel_image(~isfinite(kernel_image) | kernel_image < 0) = 0;
floor_prediction = max(mean(counts) * 1e-9, 1e-12);
history = zeros(numel(alpha), numel(record_iterations));
record_index = 1;

for iteration = 1:record_iterations(end)
    image = kernel_image + bias;
    image(~support | ~isfinite(image) | image < 0) = 0;
    predicted = max(normalization .* project(image) + background, floor_prediction);
    ratio = counts ./ predicted;
    ratio(~isfinite(ratio)) = 0;

    correction = backproject(normalization .* ratio);
    correction(~isfinite(correction) | correction < 0) = 0;
    kernel_correction = K' * correction;
    kernel_correction(~isfinite(kernel_correction) | kernel_correction < 0) = 0;
    alpha = alpha ./ kernel_sensitivity .* kernel_correction;
    alpha(~support | ~isfinite(alpha) | alpha < 0) = 0;

    kernel_image = K * alpha;
    kernel_image(~support | ~isfinite(kernel_image) | kernel_image < 0) = 0;
    bias = (reference - kernel_image) ./ (1 + lambda);
    bias(~support | ~isfinite(bias)) = 0;

    % Record the image after updating both alpha and the signed affine term.
    if iteration == record_iterations(record_index)
        image = kernel_image + bias;
        image(~support | ~isfinite(image) | image < 0) = 0;
        history(:,record_index) = image;
        record_index = record_index + 1;
    end
end
end
