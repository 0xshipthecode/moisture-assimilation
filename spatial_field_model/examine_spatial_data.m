
%
% Look at the station data with respect to current topography.
%

% Read geographic position of grid points
lon = ncread('../real-data/wrfout_d03_latlon.nc', 'XLONG');
lon = lon(:, :, 1);
lat = ncread('../real-data/wrfout_d03_latlon.nc', 'XLAT');
lat = lat(:, :, 1);

% Read in the variables
T = ncread('../real-data/wrfout_d03_T2.nc', 'T2');
Q = ncread('../real-data/wrfout_d03_Q2.nc', 'Q2');
P = ncread('../real-data/wrfout_d03_PSFC.nc', 'PSFC');

% compute moisture equilibrium field
[Ed, Ew] = equilibrium_moisture(P, Q, T);

% visualize this for all time steps
for i=1:size(Ed,3)
    subplot(121);
    imagesc(Ed(:, end:-1:1, i)');
    title(sprintf('Drying, frame %d', i));
    caxis([0, 0.55]);
    colorbar();
    subplot(122);
    imagesc(Ew(:, end:-1:1, i)');
    title(sprintf('Wetting, frame %d', i));
    caxis([0, 0.55]);
    colorbar();
    pause(0.2);
end
