function startup_fem2d()
% Add all framework subdirectories to MATLAB path.
root = fileparts(mfilename('fullpath'));
addpath(root);
addpath(fullfile(root, 'mesh'));
addpath(fullfile(root, 'elements'));
addpath(fullfile(root, 'elements', 'gauss'));
addpath(fullfile(root, 'materials'));
addpath(fullfile(root, 'assembly'));
addpath(fullfile(root, 'bc'));
addpath(fullfile(root, 'solver'));
addpath(fullfile(root, 'postproc'));
addpath(fullfile(root, 'examples'));
addpath(fullfile(root, 'tests'));
disp('fem2d framework loaded.');
end
