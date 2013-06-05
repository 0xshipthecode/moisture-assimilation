

%
%  Plots the results from test_kriging_estimators.jl
%
%  Results are csv, first line is true values, last line is no. of errors
%
%

estimators = { 'ml_bisect', 'moment_bisect', 'estimator_e_1', 'constant_diagonal' };
Nest = length(estimators);
NNobs = 3;

fprintf('| Estimator  | Nobs  | mean $\\hat{\\sigma}^2$  | est. bias $\\hat{\\sigma}^2 - \\sigma^2$ | std $\\hat{\\sigma}^2$  | coeff of var. $\\hat{\\sigma}^2$  | sum(I($\\hat{\\sigma}^2$ < 0)) |\n');

figure;
max_val = 0;
min_val = 0;
ndx = 1;    
for i = 1:Nest
    est = estimators{i};

    fprintf('| %s     |   |   |   |   |   |\n', est);
    for nobs = [ 20, 50, 100 ]
        
        fname = sprintf('test-%s-%d.log', est, nobs);
        
        t = load(fname);
        t_s2 = t(2:end-1, 3);
        
        subplot(Nest, NNobs, ndx);
        hist(t_s2, 45);
        min_val = min([min_val; t_s2]);
        max_val = max([max_val; t_s2]);
        a = axis();
        hold on
        plot([t(1,3), t(1,3)], [0, a(4)], 'r-', 'linewidth', 2);
        hold off
        
        title(sprintf('%s / %d obs', regexprep(est, '_', '-'), nobs));
        ndx = ndx + 1;
        
        fprintf('      |  %d | %g | %g | %g | %g | %d |\n', ...
            nobs, mean(t_s2), mean(t_s2) - t(1,3), std(t_s2), std(t_s2)/mean(t_s2), t(end,1));
        
    end
    
end


ndx = 1;
for i=1:Nest
    for j=1:NNobs
        subplot(Nest, NNobs, ndx);
        a = axis();
        axis([min_val, max_val, a(3), a(4)]);
        ndx = ndx + 1;
    end
end
        

% load test1.log
% t1 = test1;
% clear test1
% 
% load test2.log
% t2 = test2;
% clear test2
% 
% load test3.log
% t3 = test3;
% clear test3
% 
% titles = { '$\beta_1$', '$\beta_2$', '$\sigma^2_\eta$' };
% 
% 
% figure;
% for i=1:3
%     subplot(330+i);
%     hist(t1(2:end-1, i), 45);
%     a = axis();
%     hold on
%     plot([t1(1,i), t1(1,i)], [0, a(4)], 'r-', 'linewidth', 2);
%     hold off
%     title(sprintf('Estimator 1 [%s]', titles{i}), 'fontsize', 18, 'interpreter', 'latex');
% end
% 
% for i=1:3
%     subplot(333+i);
%     hist(t2(2:end-1, i), 45);
%     a = axis();
%     hold on
%     plot([t2(1,i), t2(1,i)], [0, a(4)], 'r-', 'linewidth', 2);
%     hold off
%     title(sprintf('Estimator 2 [%s]', titles{i}), 'fontsize', 18, 'interpreter', 'latex');
% end
% 
% 
% for i=1:3
%     subplot(336+i);
%     hist(t3(2:end-1, i), 45);
%     a = axis();
%     hold on
%     plot([t3(1,i), t3(1,i)], [0, a(4)], 'r-', 'linewidth', 2);
%     hold off
%     title(sprintf('Estimator 3 [%s]', titles{i}), 'fontsize', 18, 'interpreter', 'latex');
% end
% 
% % run through plots and adjust axes
% for i=1:3
%     subplot(330+i);
%     a1 = axis();
%     subplot(333+i);
%     a2 = axis();
%     subplot(336+i);
%     a3 = axis();
%     
%     a = [ min([a1(1), a2(1), a3(1)]), max([a1(2), a2(2), a3(2)]), ...
%           min([a1(3), a2(3), a3(3)]), max([a1(4), a2(4), a3(4)]) ];
%       
%     subplot(330+i);
%     axis(a);
%     subplot(333+i);
%     axis(a);
%     subplot(336+i);
%     axis(a);
% end
% 
% % print some summary statistics
% t1_s2 = t1(2:end-1, 3);
% t2_s2 = t2(2:end-1, 3);
% t3_s2 = t3(2:end-1, 3);
% 
% fprintf('Estimator 1: mean %g   std %g   coeff_of_var %g\n', mean(t1_s2), std(t1_s2), std(t1_s2)/mean(t1_s2));
% fprintf('Estimator 2: mean %g   std %g   coeff_of_var %g\n', mean(t2_s2), std(t2_s2), std(t2_s2)/mean(t2_s2));
% fprintf('Estimator 3: mean %g   std %g   coeff_of_var %g\n', mean(t3_s2), std(t3_s2), std(t3_s2)/mean(t3_s2));
