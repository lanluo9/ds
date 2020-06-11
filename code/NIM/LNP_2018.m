% 2018-09-26-0/ data007-nyj-map/ white noise; data005 (195-4000s) & data006 flow 
clear
clc
close

datapath = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/2018-09-26-0/data007-nyj-map/data007-nyj-map';
% datapath = '/Volumes/dusom_fieldlab/All_Staff/lab/Experiments/Array/Analysis/2018-09-26-0/data007-nyj-map/data007-nyj-map';
% datapath = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2018-09-26-0/data007-nyj-map/data007-nyj-map';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);

%%
% movie_path = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/ds/code/NIM/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
% % movie_path = '/Users/circuit/Documents/MATLAB/matlab/private/Lan/ds/code/NIM/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';
% % movie_path = '/Users/katthanya/Documents/MATLAB/matlab/private/Lan/code/NIM/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml';

% % mvi = load_movie(movie_path, datarun.triggers);
% [mov, ~,~, dur, refresh] = get_movie_LL(movie_path, datarun.triggers, 216000); 
% mov = squeeze(mov(:,:,1,:)); % no need for color dimension 
% [NY, NX, NFRAMES] = size(mov);
% 
% save('mov_20180926_007.mat', 'mov', '-v7.3') % force save >2GB .mat

%% 
% cell_index_array = [4:13, find(datarun.cell_ids == 1578), find(datarun.cell_ids == 4756)];
% for counter = length(cell_index_array) : -1 : 1

cell_index = 3; %cell_index_array(counter);
cell_id = datarun.cell_ids(cell_index)
spikes = datarun.spikes{cell_index};

% % reduce Xstim size in time
% re_time = 216000 / 60 * 10; % take first 10 mins
% mov_resize = mov(:, :, 1:re_time);
% mov = reshape(mov_resize, [NY*NX, re_time])';

% reduce Xstim size in space
% cd /Users/katthanya/Documents/MATLAB/matlab/private/Lan/code/NIM
cd D:\RRR\Grad\Rotation\GF_lab\lab_Mac\ds\code\NIM
load mov_20180926_007.mat
[NY, NX, NFRAMES] = size(mov);

RF = datarun.vision.sta_fits{cell_index,1};
RF_center_x = RF.mean(1);
RF_center_y = NY - RF.mean(2);
re_size = 13/2 - 1; % reduced from 15 bc PC has low memory

x_upper = ceil(RF_center_x + re_size);
x_lower = floor(RF_center_x - re_size);
y_upper = ceil(RF_center_y + re_size);
y_lower = floor(RF_center_y - re_size);

mov_resize = mov(y_lower:y_upper, x_lower:x_upper, :);
NX = size(mov_resize,1); NY = NX; % redefine size after resizing
mov = reshape(mov_resize, [NY*NX, NFRAMES])';
mov(mov < 0.5) = -0.48; % taken from xml contrast value
mov(mov > 0.5) = 0.48;

%% 

up_samp_fac = 1; % 1; A little nicer resolution at 2, but indeed runs longer
tent_basis_spacing = 1;
nLags = 25*up_samp_fac; % 25; 
dt = 16.5975 ./ 1000;
binsize = dt/up_samp_fac;
NT = NFRAMES*up_samp_fac;

% Init stim parameters structure
params_stim = NIM.create_stim_params([nLags NY, NX], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);

% generate stimulus matrix
Xstim = NIM.create_time_embedding(mov, params_stim); % beware of memory
% bin spike times 
Robs = NIM.Spks2Robs(spikes, binsize, NT );

% For generating nested-cross validation indices (this is default 5-fold)
[Ui, XVi] = NIM.generate_XVfolds( NT ); % XVi array size = NT/5. XVi is test set, while Ui is training set


%% round 1 fit
disp('Fitting GLM: 1 excitatory filter')
mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)
NL_types = {'lin'}; % define subunit as linear 

% New object-based definition of model (note regularization can be set here)
fit0 = NIM( params_stim, NL_types, mod_signs, 'd2xt', 100 ); 
% also using spatial regularization too -- can do seprately or as one

% Tweak optimization parameters to doesn't creep around the minimum forever (necessary for 
% some reason with this data)
optim_params.maxIter = 400;

% Fit stimulus filters
fit0 = fit0.fit_filters( Robs, Xstim, Ui, 'optim_params', optim_params );

% find optimal regularization using nested cross-val (note this can take a while)
% fit_test = fit0.reg_path( Robs, Xstim, Ui, XVi, 'lambdaID', 'd2xt', 'optim_params', optim_params );
% Subunit 1: Best reg = 600.00 for reduction in space mov; 400 in time; 200
% for spatial RF
% this finds optimal value of XX

% look at filter
k0 = reshape(fit0.subunits(1).filtK, [nLags, NY * NX]);
[mval,bestlag] = max(max(abs(k0)')); bestlag 


figure(2); clf;
[U, S, V] = svd(k0);
subplot(1,2,1)
plot(U(:,1))
axis tight
axis square
subplot(1,2,2); colormap gray
imagesc(reshape(V(:,1), NY, NX))
axis square
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['k0_filter-' num2str(cell_id) '.png'])
close

% %% round 2 fit NIM
mod_signs = [1 -1]; % both inputs are excitatory. (+1 for excitatory, -1 for suppressive)
NL_types = {'rectlin','rectlin'}; % name-change

% Try spike history term
fitS = fit0.init_spkhist(16,'doubling_time',4, 'negcon'); % I don't believe in positive spk-history
fitS = fitS.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );

% Always best to seed inhibition with different possibilities: random will often get stuck in local minima
fit1 = fit0;
fit1.subunits(1).NLtype='rectlin'; % note new name (from threshlin)
fit1.subunits(2) = fit1.subunits(1); % make second subunit like the first
fit1.subunits(2).weight = -1; % make suppressive
kshift = NIM.shift_mat_zpad(k0,4,1); % shift filter to be delayed
fit1.subunits(2).filtK = kshift(:); % flatten

fit1 = fit1.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
%fit1 = fit1.reg_path( Robs, Xstim, Ui, XVi, 'lambdaID', 'd2xt', 'fit_offsets', 1, 'optim_params', optim_params );

% Fit with spike-history term
fit1S = fit1.init_spkhist(16,'doubling_time',4, 'negcon');
fit1S = fit1S.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
% spike history term helps
fit1S.display_model('Xstims', Xstim )
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit1S-' num2str(cell_id) '.png'])
close
% plot by itself
% plot(fit1S.spk_hist.bin_edges(2:end), fit1S.spk_hist.coefs)
% close


% new model display
fit1.display_model('Xstims',Xstim)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit1-' num2str(cell_id) '.png'])
close

% can also plot RFs separately
k1 = reshape(fit1.subunits(1).filtK, [nLags, NY*NX]);
k2 = reshape(fit1.subunits(2).filtK, [nLags, NY*NX]);
[~,best_pixel] = max(max(abs(k1)));
[mk1,blat1] = max(max(abs(k1')));
[mk2,blat2] = max(max(abs(k2')));

figure
subplot(1,3,1); hold on
plot(k1(:,best_pixel)/std(k1(:,best_pixel)),'b')
plot(k2(:,best_pixel)/std(k2(:,best_pixel)),'r')
axis tight
% Note that there is really only a delay between the filters of ~4 ms, which means that they will have a hard
% time explaining differences in timing at the resolution of the PSTH (where event width >>4 ms)
subplot(1,3,2); colormap gray
imagesc(reshape(k1(blat1,:),NY,NX),[-1 1]*mk1)
axis square
subplot(1,3,3); colormap gray
imagesc(reshape(k2(blat2,:),NY,NX), [-1 1]*mk2)
axis square
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['subunits-' num2str(cell_id) '.png'])
close

fit0nl = fit0.fit_spkNL( Robs, Xstim, Ui );
fit1nl = fit1.fit_spkNL( Robs, Xstim, Ui );

% fit upstream nonlinearity?
fit2 = fit1.init_nonpar_NLs( Xstim );
fit2 = fit2.fit_upstreamNLs( Robs, Xstim, Ui );
% this actually ends up looking rectlin, so not worth it
fit2.display_model('Xstims',Xstim)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit2-' num2str(cell_id) '.png'])
close

% internal LLx (cross-validated likelihood, although was used for meta-params)
[LLs,~,~,fp] = fit0.eval_model(Robs, Xstim, XVi );
LLs(2) = fitS.eval_model(Robs, Xstim, XVi );
LLs(3) = fit1.eval_model(Robs, Xstim, XVi );
LLs(4) = fit1S.eval_model(Robs, Xstim, XVi );
LLs(5) = fit2.eval_model(Robs, Xstim, XVi ); % fit2 (w nonlin) is not better than fit1
LLfit = LLs - fp.nullLL

% save('before_refit_63_13.mat', '-regexp', 'fit*')

%%
% load before-refit-63.mat 

[~,pred_rate0,~,~] = fit0.eval_model(Robs, Xstim, XVi );
[~,pred_rateS,~,~] = fitS.eval_model(Robs, Xstim, XVi );
[~,pred_rate1,~,~] = fit1.eval_model(Robs, Xstim, XVi );
[~,pred_rate1S,~,~] = fit1S.eval_model(Robs, Xstim, XVi );
[~,pred_rate2,~,~] = fit2.eval_model(Robs, Xstim, XVi );
pred_rates = [pred_rate0, pred_rateS, pred_rate1, pred_rate1S, pred_rate2];

Robs_test = Robs(XVi); % pred_rate unit is spike per bin, not per sec. pred_rate array correspond to test set
Robs_test_shuffled = Robs_test(randperm(length(Robs_test))); % shuffle spk to visualize performance

% %% plot Robs real spk vs pred_rate
% color = prism(size(pred_rates,2));
% bgn = 42900;
% fin = 43000;
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% for i = 1:size(pred_rates,2)
%     pred_rate_now = pred_rates(:, i);
%     plot_pred = plot(pred_rate_now(bgn:fin), 'Color', color(i,:), 'LineWidth', 1);
%     plot_pred.Color(4) = 0.7; % alpha transparency
%     hold on
% end
% plot_spk = plot(Robs_test(bgn:fin), 'k', 'LineWidth', 2);
% plot_spk.Color(4) = 0.4; 
% legend({'0','S','1', '1S', '2', 'robs'}, 'Location','northeast')
% legend('boxoff')
% xlim([0-5, (fin-bgn)+5])
% saveas(gcf, ['perf-spk-' num2str(cell_id) '.png'])

%% compare shuffled Robs vs real Robs vs pred_rate
color = prism(size(pred_rates,2));
bgn = 200;
fin = 400;

subplot(2,1,1)
for i = 1:size(pred_rates,2)
    pred_rate_now = pred_rates(:, i);
    plot_pred = plot(bgn:fin, pred_rate_now(bgn:fin), 'Color', color(i,:), 'LineWidth', 1);
    plot_pred.Color(4) = 0.7; % alpha transparency
    hold on
end
plot_spk = plot(bgn:fin, Robs_test(bgn:fin), 'k', 'LineWidth', 2);
plot_spk.Color(4) = 0.4; 
legend({'0','S','1', '1S', '2', 'robs'}, 'Location','northeast')
legend('boxoff')

subplot(2,1,2)
for i = 1:size(pred_rates,2)
    pred_rate_now = pred_rates(:, i);
    plot_pred = plot(bgn:fin, pred_rate_now(bgn:fin), 'Color', color(i,:), 'LineWidth', 1);
    plot_pred.Color(4) = 0.7; % alpha transparency
    hold on
end
plot_spk = plot(bgn:fin, Robs_test_shuffled(bgn:fin), 'k', 'LineWidth', 2);
plot_spk.Color(4) = 0.4; 

set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['perf-spk-cp-' num2str(cell_id) '.png'])


%% test R2
for i = 1:5
    R2(i) = 1 - mean( (Robs_test - pred_rates(:,i)) .^2) / var(Robs_test);
    R2_shuffled(i) = 1 - mean( (Robs_test_shuffled - pred_rates(:,i)) .^2) / var(Robs_test_shuffled);
end
R2
R2_shuffled
sum(R2>0) / size(R2,1) / size(R2,2) 
sum(R2_shuffled>0) / size(R2_shuffled,1) / size(R2_shuffled,2) 

% save('before_refit_63_13.mat', 'LLfit', 'R2', '-append')

%% plot PSTH instead of Robs real spk % ofc this should not be done
% 
% rep_len = 50; % pretend there are reps of white noise. unit = frames. 1s = 60 frames
% Robs_test_rep = reshape(Robs_test, [rep_len, length(Robs_test)/rep_len]);
% Robs_test_rep = sum(Robs_test_rep, 1) ./ size(Robs_test_rep, 2); 
% 
% color = prism(size(pred_rates,2));
% bgn = 50;
% fin = 150;
% 
% for i = 1:size(pred_rates,2)
%     pred_rate_now = pred_rates(:, i);
%     pred_rate_now = reshape(pred_rate_now, [rep_len, length(pred_rate_now)/rep_len]);
%     pred_rate_now = sum(pred_rate_now, 1) ./ size(pred_rate_now, 2);    
%     
%     plot_pred = plot(pred_rate_now(bgn:fin), 'Color', color(i,:), 'LineWidth', 1);
%     plot_pred.Color(4) = 0.7; % alpha transparency
%     hold on
% end
% plot_spk = plot(Robs_test_rep(bgn:fin), 'k', 'LineWidth', 2);
% plot_spk.Color(4) = 0.4; 
% 
% legend({'0','S','1', '1S', '2', 'robs'}, 'Location','northeast')
% legend('boxoff')
% xlim([0-2, (fin-bgn)+5])
% set(gcf, 'Position', get(0, 'Screensize'));
% saveas(gcf, ['perf-PSTH-' num2str(cell_id) '.png'])
% 
% % save(['perf-' num2str(cell_id) '.mat'], '-v7.3')
% % save('perf_63_noXstim.mat', '-regexp', '^((?!Xstim).)*$') % not serializable?
% disp(['finished eval_model for cell' num2str(cell_id)])
% 
% % toc
% % end