%
% this script automatically sets up the paths (if needed) to
% use the moisture assimilationt toolkit
%

if(~exist('ukf_select_sigma_points.m','file'))
  addpath(genpath('.'));
end