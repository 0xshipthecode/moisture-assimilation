%% load all the examine_mesowest_stations.py output
d = load('distance_matrix.txt');
wd = load('wdistance_matrix.txt');
e = load('elevations.txt');
f = load('fm10_data.txt');
T = load('T_data.txt');
wf = load('wfm10_data.txt');
wT = load('wT_data.txt');

%% extract cov or corr matrices
rel = @corrcoef;
rel_str = 'cc';

% rel = @cov;
% rel_str = 'cov';

%sel = [  3, 5, 7, 8, 9, 11, 15, 16, 22, 25, 28, 30, 32, 33];

c = rel(f);
wc = rel(wf);
cT = rel(T);
wcT = rel(wT);


%% extract and vectorize upper triangular parts
utri = triu(true(38),1);
dt = d(utri);
wdt = wd(utri);
et = e(utri);
ct = c(utri);
wct = wc(utri);
cTt = cT(utri);
wcTt = wcT(utri);

% show scatter plots
figure;
scatter(ct,wct);
title('Scatter of fm10 covariance vs. nearest grid point fm10 covar');
saveas(gcf, 'fm10_station_vs_model_covar.png', 'png');

figure;
scatter(cTt, wcTt);
title('Scatter of T2 covar vs. nearest grid point T2 covar');
saveas(gcf, 'T2_station_vs_model_covar.png', 'png');

% show rel vs. distance plots
figure;
plot(dt, ct, 'ro', wdt, wct, 'go');
title(sprintf('%s vs. distance for stations and for model fm10', rel_str));
legend('station', 'wrf');
saveas(gcf, sprintf('%s_vs_dist_fm10.png', rel_str), 'png');

figure;
plot(dt, cTt, 'ro', wdt, wcTt, 'go');
title(sprintf('%s vs. distance for stations and for model T2', rel_str));
legend('station', 'wrf');
saveas(gcf, sprintf('%s_vs_dist_T2.png', rel_str), 'png');

% Show some correlations
fprintf('Correlation/Covariances between model and experimental variables\n')
fprintf('Temp: pearson  %g  spearman %g\n', corr(cTt,wcTt), corr(cTt,wcTt,'type','spearman'));
fprintf('fm10: pearson  %g  spearman %g\n', corr(ct,wct), corr(ct,wct,'type','spearman'));
fprintf('distance: pearson  %g  spearman %g\n', corr(dt,wdt), corr(dt,wdt,'type','spearman'));

