%DEMO_ZUBAL_AFFINEKEM Reconstruct and display Zubal frame 12.
% Set FESSLER_IRT_ROOT and KER_ROOT before running this script.

demo_results = run_demo(fileparts(fileparts(mfilename('fullpath'))));

function results = run_demo(root_dir)
previous_path = path;
restore_path = onCleanup(@() path(previous_path)); %#ok<NASGU>
addpath(root_dir);

results = run_zubal_affinekem(Frames=12, Iterations=30, PriorIterations=20, ...
    OutputDir=fullfile(root_dir, 'results', 'demo'));
final_metrics = results.metrics(end,:);
disp(final_metrics);

panels = {results.truth(:,:,1), results.reference(:,:,1), results.images(:,:,end,1)};
panel_titles = {'Ground truth', 'Composite reference', 'AffineKEM'};
intensity_limits = [0, max(panels{1}(:))];

fig = figure('Name', 'Zubal AffineKEM demo', 'Color', 'w');
layout = tiledlayout(fig, 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('Zubal frame %d | %d iterations | NMSE %.2f dB', ...
    final_metrics.Frame, final_metrics.Iteration, final_metrics.NMSE_dB));
colormap(fig, gray);

for index = 1:numel(panels)
    ax = nexttile(layout);
    imagesc(ax, panels{index}, intensity_limits);
    axis(ax, 'image', 'off');
    title(ax, panel_titles{index});
end

scale_bar = colorbar(ax);
scale_bar.Layout.Tile = 'east';
end
