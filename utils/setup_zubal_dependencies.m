function setup_zubal_dependencies(irt_path, ker_path)
%SETUP_ZUBAL_DEPENDENCIES Load the external IRT and KER utilities.
assert(strlength(irt_path) > 0 && isfile(fullfile(irt_path,'setup.m')), ...
    'Set IRTPath or FESSLER_IRT_ROOT to the IRT directory containing setup.m.');
assert(strlength(ker_path) > 0 && isfile(fullfile(ker_path,'utils','buildSparseK.m')), ...
    'Set KERPath or KER_ROOT to the KER v0.2 directory containing utils/buildSparseK.m.');
% Resolve the root before IRT setup changes its working directory.
[ok, info] = fileattrib(char(irt_path));
assert(ok, 'Could not resolve the IRT directory.');
irtdir = info.Name; %#ok<NASGU>
run(fullfile(irtdir, 'setup.m'));
addpath(fullfile(ker_path, 'utils'));
end
