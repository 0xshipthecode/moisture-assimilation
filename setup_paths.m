%
% this script automatically sets up the paths (if needed) to
% use the moisture assimilationt toolkit
%

if(~exist('ukf_select_sigma_points.m','file'))
  cp = pwd();
  addpath(fullfile(cp, 'extended_model'));
  addpath(fullfile(cp, 'state_model'));
  addpath(fullfile(cp, 'tests'));
  addpath(fullfile(cp, 'ukf'));
end