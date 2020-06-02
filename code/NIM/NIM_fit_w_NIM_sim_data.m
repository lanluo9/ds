%%
clc
clear
close all

% load perf-1578.mat
load demo_adapt.mat
%% 
[~,pred_rate0,~,~] = fit0.eval_model(Robs, Xstim, XVi );
[~,pred_rateS,~,~] = fitS.eval_model(Robs, Xstim, XVi );
[~,pred_rate1,~,~] = fit1.eval_model(Robs, Xstim, XVi );
[~,pred_rate1S,~,~] = fit1S.eval_model(Robs, Xstim, XVi );
[~,pred_rate2,~,~] = fit2.eval_model(Robs, Xstim, XVi );

%%
pred_rates = [pred_rate0, pred_rateS, pred_rate1, pred_rate1S, pred_rate2];

Robs_test = Robs(XVi);
% Robs_5 = reshape(Robs, [size(pred_rates,2), size(pred_rates,1)]);
% dt = 16.5975 ./ 1000;
% Robs_5 = sum(Robs_5, 1) ./ 5; % Robs_5 is calculated for every 5 frames

% edges = linspace(floor(spikes(1)), ceil(spikes(end)), size(pred_rates,1) + 1 );
% [spk_binned, ~] = histcounts(spikes, edges);
% spk_per_bin = spk_binned ./ 5; % pred_rate is spike per bin, not per sec
% % spk_per_sec = spk_binned ./ (binsize * 5); 

%%
color = prism(size(pred_rates,2));
bgn = 300;
fin = 500;

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:size(pred_rates,2)
    pred_rate_now = pred_rates(:, i);
    plot_pred = plot(pred_rate_now(bgn:fin), 'Color', color(i,:), 'LineWidth', 1);
    plot_pred.Color(4) = 0.7; % alpha transparency
    hold on
end
plot_spk = plot(Robs_test(bgn:fin), 'k', 'LineWidth', 2);
plot_spk.Color(4) = 0.4; 
legend({'0','S','1', '1S', '2', 'robs'}, 'Location','northeast')
legend('boxoff')
xlim([0-5, (fin-bgn)+5])

%% test R2
% R2_part = 1 - mean( (spk_per_bin(bgn:fin)' - pred_rates(bgn:fin,4)) .^2) / var(spk_per_bin(bgn:fin))
for i = 1:5
    R2(i) = 1 - mean( (Robs_test - pred_rates(:,i)) .^2) / var(Robs_test);
end
R2
% plot(R2)
sum(R2>0) / size(R2,1) / size(R2,2) % model display seems reasonable though, should be due to buggy eval?

%% plot PSTH instead of Robs
% 1s = 60 frames

rep_len = 50; % pretend there are reps of white noise. unit = frames
Robs_test_rep = reshape(Robs_test, [rep_len, length(Robs_test)/rep_len]);
Robs_test_rep = sum(Robs_test_rep, 1) ./ size(Robs_test_rep, 2); 

color = prism(size(pred_rates,2));
bgn = 50;
fin = 150;

figure('units','normalized','outerposition',[0 0 1 1])
for i = 1:size(pred_rates,2)
    pred_rate_now = pred_rates(:, i);
    pred_rate_now = reshape(pred_rate_now, [rep_len, length(pred_rate_now)/rep_len]);
    pred_rate_now = sum(pred_rate_now, 1) ./ size(pred_rate_now, 2);    
    
    plot_pred = plot(pred_rate_now(bgn:fin), 'Color', color(i,:), 'LineWidth', 1);
    plot_pred.Color(4) = 0.7; % alpha transparency
    hold on
end
plot_spk = plot(Robs_test_rep(bgn:fin), 'k', 'LineWidth', 2);
plot_spk.Color(4) = 0.4; 
legend({'0','S','1', '1S', '2', 'robs'}, 'Location','northeast')
legend('boxoff')
xlim([0-2, (fin-bgn)+5])
% saveas(gcf, ['perf-' num2str(cell_id) '.png'])



%% generate fake spike train from model
% r_array = poissrnd(20,[1000 1]);
% mean(r_array)
% 
% cell_id = 63;
% save(['perf-test-' num2str(cell_id) '.mat'], 'fit0', '-v7.3')



