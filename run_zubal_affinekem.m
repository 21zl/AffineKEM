function results = run_zubal_affinekem(options)
%RUN_ZUBAL_AFFINEKEM Simulate dynamic Zubal PET and reconstruct selected frames.
arguments
    options.IRTPath (1,1) string = string(getenv('FESSLER_IRT_ROOT'))
    options.KERPath (1,1) string = string(getenv('KER_ROOT'))
    options.DataFile (1,1) string = ""
    options.OutputDir (1,1) string = "results"
    options.Frames (1,:) double {mustBeInteger,mustBeInRange(options.Frames,1,24)} = [2 12 24]
    options.Iterations (1,1) double {mustBeInteger,mustBePositive} = 100
    options.PriorIterations (1,1) double {mustBeInteger,mustBePositive} = 100
    options.RecordEvery (1,1) double {mustBeInteger,mustBePositive} = 10
    options.TotalCounts (1,1) double {mustBePositive,mustBeFinite} = 8e6
    options.Seed (1,1) double {mustBeInteger,mustBeInRange(options.Seed,0,4294967295)} = 2026080900
    options.Lambda (1,1) double {mustBeNonnegative,mustBeFinite} = 0.8
    options.ReferenceSigma (1,1) double {mustBeNonnegative,mustBeFinite} = 0.5
end

frames = options.Frames;
assert(~isempty(frames) && numel(unique(frames)) == numel(frames), ...
    'Frames must contain distinct frame numbers from 1 to 24.');

root = fileparts(mfilename('fullpath'));
if strlength(options.DataFile) == 0
    options.DataFile = fullfile(root, 'data', 'zubal_data.mat');
    if ~isfile(options.DataFile)
        options.DataFile = fullfile(options.KERPath, 'demo', 'zubal_data.mat');
    end
end

previous_rng = rng;
restore_rng = onCleanup(@() rng(previous_rng)); %#ok<NASGU>
previous_path = path;
restore_path = onCleanup(@() path(previous_path)); %#ok<NASGU>
addpath(fullfile(root, 'functions'), fullfile(root, 'utils'));
setup_zubal_dependencies(options.IRTPath, options.KERPath);
rng(options.Seed, 'twister');

[counts, normalization, background, truth, durations, image_size, project, backproject] = ...
    simulate_zubal_pet(options.DataFile, options.TotalCounts);
[K, composite_reference, group_for_frame, group_duration] = build_affinekem_priors( ...
    project, backproject, counts, normalization, background, durations, image_size, options.PriorIterations);

iterations = unique([options.RecordEvery:options.RecordEvery:options.Iterations, options.Iterations]);
n_frames = numel(frames);
results.images = zeros([image_size, numel(iterations), n_frames]);
results.reference = zeros([image_size, n_frames]);
results.truth = reshape(truth(:,frames), [image_size, n_frames]);
mse = zeros(numel(iterations), n_frames);
nmse_db = zeros(size(mse));

for index = 1:n_frames
    frame = frames(index);
    fprintf('Reconstructing frame %d (%d/%d).\n', frame, index, n_frames);
    group = group_for_frame(frame);
    ref = composite_reference(:,group) / group_duration(group) * durations(frame);
    ref = reshape(ref, image_size);
    if options.ReferenceSigma > 0
        ref = imgaussfilt(ref, options.ReferenceSigma, 'Padding', 'replicate');
    end
    results.reference(:,:,index) = ref;
    history = affinekem_reconstruct(project, backproject, K, counts(:,frame), ...
        normalization, background(:,frame), ref(:), options.Lambda, iterations);
    results.images(:,:,:,index) = reshape(history, [image_size, numel(iterations)]);
    squared_error = (history - truth(:,frame)).^2;
    mse(:,index) = mean(squared_error, 1).';
    nmse_db(:,index) = 10*log10(sum(squared_error, 1).' / sum(truth(:,frame).^2));
end

[iteration_grid, frame_grid] = ndgrid(iterations, frames);
results.metrics = table(frame_grid(:), iteration_grid(:), mse(:), nmse_db(:), ...
    'VariableNames', {'Frame','Iteration','MSE','NMSE_dB'});
results.frames = frames;
results.iterations = iterations;
results.frame_duration_minutes = durations(frames);
results.settings = rmfield(options, {'IRTPath','KERPath','DataFile','OutputDir'});
results.settings.UsePSF = false;

if strlength(options.OutputDir) > 0
    if ~isfolder(options.OutputDir)
        mkdir(options.OutputDir);
    end
    save(fullfile(options.OutputDir, 'zubal_affinekem.mat'), 'results');
    writetable(results.metrics, fullfile(options.OutputDir, 'zubal_affinekem_metrics.csv'));
    fprintf('Saved AffineKEM results to %s.\n', options.OutputDir);
end
end
