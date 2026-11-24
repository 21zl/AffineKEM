# AffineKEM

PET image reconstruction using affine kernel methods.

This MATLAB example generates a 24-frame dynamic PET acquisition from the 2-D Zubal head phantom, adds Poisson noise, and reconstructs selected frames with AffineKEM. The default example reconstructs frames 2, 12, and 24.

## Repository layout

- [`run_zubal_affinekem.m`](run_zubal_affinekem.m) runs the full example and writes the results.
- [`demo/demo_zubal_affinekem.m`](demo/demo_zubal_affinekem.m) runs and displays a shorter reconstruction of frame 12.
- [`functions/simulate_zubal_pet.m`](functions/simulate_zubal_pet.m) generates the dynamic ground truth and the noisy PET projections.
- [`functions/build_affinekem_priors.m`](functions/build_affinekem_priors.m) builds the sparse kernel and composite reference images.
- [`functions/affinekem_reconstruct.m`](functions/affinekem_reconstruct.m) reconstructs one frame.
- [`utils/setup_zubal_dependencies.m`](utils/setup_zubal_dependencies.m) configures the external toolbox paths. The entry point adds the local folders automatically.

## Requirements

- MATLAB R2022a with Image Processing Toolbox and Statistics and Machine Learning Toolbox.
- Jeffrey A. Fessler's Michigan Image Reconstruction Toolbox (IRT/MIRT). Install a version with a working `Gtomo2_strip` projector and the MEX binaries for your MATLAB platform. [MIRT download and setup](https://web.eecs.umich.edu/~fessler/code/index.html)
- Guobao Wang's KER v0.2. Install KER separately for its sparse kernel construction utilities in `utils/`. [Wang Lab software downloads](https://wanglab.faculty.ucdavis.edu/code)

The Zubal phantom file ships with the KER v0.2 demo. The entry point looks for a local copy at `data/zubal_data.mat` and then for `demo/zubal_data.mat` under `KER_ROOT`. To use a file stored elsewhere, pass `DataFile="/path/to/zubal_data.mat"`. The MAT file must contain the `111 × 111` label and attenuation maps (`model` and `u`), the `24 × 6` activity curves (`TAC`), and the `24 × 2` frame start/end times in seconds (`scant`).

## Run

Open MATLAB in this directory and set the paths to the two installations.

```matlab
setenv('FESSLER_IRT_ROOT', '/path/to/irt');       % Contains setup.m
setenv('KER_ROOT', '/path/to/KER_v0.2');        % Contains utils/
run('demo/demo_zubal_affinekem.m');
```

The demo uses 20 composite image iterations and 30 AffineKEM iterations. It displays the ground truth, composite reference, and reconstruction on a shared intensity scale, and writes the MAT and CSV files under `results/demo/`.

For the full three-frame experiment, run:

```matlab
results = run_zubal_affinekem();
```

Run a quick check with a single frame and two iterations.

```matlab
quick = run_zubal_affinekem(Frames=12, Iterations=2, PriorIterations=2, ...
    OutputDir="");
```

The default run writes `results/zubal_affinekem.mat` and `results/zubal_affinekem_metrics.csv` under the current directory. The MAT file contains the reconstructed images at the recorded iterations, the corresponding ground truth and affine references, and the run settings. Images are indexed as `results.images(:,:,iterationIndex,frameIndex)`. Use `results.iterations` and `results.frames` to map the indices to iteration and frame numbers. The CSV reports MSE and normalized MSE in dB over the full image. NMSE in dB is defined as `10*log10(sum((image-truth).^2)/sum(truth.^2))` over the pixels.

To view frame 12 at the final recorded iteration of the default run:

```matlab
imagesc(results.images(:,:,end,2));
axis image off;
colormap gray;
colorbar;
```

## Simulation and reconstruction settings

| Setting | Default |
| --- | --- |
| Acquisition | 2-D parallel-beam strip projector with attenuation, uniform background at 20% of the mean attenuated signal per frame, and Poisson noise |
| Image and sinogram | 111 × 111 pixels and 249 bins × 210 angles |
| Expected total counts | 8 million across 24 frames, including background |
| Poisson seed | 2026080900, Twister |
| Composite groups | Frames 1–16, 17–20, and 21–24 |
| Composite image estimation | 100 ML-EM iterations per image |
| Kernel | 48 nearest neighbors in composite feature space, radial width 1 |
| Affine reference | Duration-scaled composite image, Gaussian smoothing at 0.5 pixels with replicate padding |
| Affine regularization | `Lambda=0.8` |
| AffineKEM iterations | 100 per frame, recorded every 10 iterations and at the final iteration |
| System PSF | Disabled |

Composite ML-EM with summed normalization supplies the smoothed kernel features. A second composite estimate uses mean normalization to form the fixed reference, which is scaled to the target frame duration before smoothing. AffineKEM initializes the signed additive term to zero and alternates a multiplicative coefficient update with the residual update `b = (reference - K*alpha)/(1 + Lambda)`. The saved iteration image uses the updated coefficient and additive term.
