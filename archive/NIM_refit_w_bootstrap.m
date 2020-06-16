%%
clc
clear

load demo_adapt.mat Robs Xstim XVi fit0 fitS fit1 fit1S fit2
% save(['before-refit-' num2str(cell_id) '.mat'], 'fit0', 'fit0nl', 'fit1', 'fit1nl', ...
%     'fit1S', 'fit2', 'fitS', '-v7.3')

[~,pred_rate0,~,~] = fit0.eval_model(Robs, Xstim, XVi );
[~,pred_rateS,~,~] = fitS.eval_model(Robs, Xstim, XVi );
[~,pred_rate1,~,~] = fit1.eval_model(Robs, Xstim, XVi );
[~,pred_rate1S,~,~] = fit1S.eval_model(Robs, Xstim, XVi );
[~,pred_rate2,~,~] = fit2.eval_model(Robs, Xstim, XVi );

%% generate fake binned spike from model
pred_rates = [pred_rate0, pred_rateS, pred_rate1, pred_rate1S, pred_rate2];

spk_bootstrap = round(pred_rates); % bootstrapping when generating fake data is equivalent to round(rate)
spk_rnd = pi * ones(size(pred_rates));
for i = 1 : size(pred_rates, 2)
    for j = 1 : size(pred_rates, 1)
        spk_rnd(j,i) = poissrnd(pred_rates(j,i));
    end
end
spk_rnd(spk_rnd < 0) = 0; % spk number cannot be negative

%% re-fit model w fake data

clear Robs Xstim XVi 
clc

datapath = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2018-09-26-0/data007-nyj-map/data007-nyj-map';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);

cell_index = 3; % cell_id = 63
cell_id = datarun.cell_ids(cell_index)

%% 

cd /Users/katthanya/Documents/MATLAB/matlab/private/Lan/code/NIM
load mov_20180926_007.mat
[NY, NX, NFRAMES] = size(mov);

up_samp_fac = 1; % 1; A little nicer resolution at 2, but indeed runs longer
tent_basis_spacing = 1;
nLags = 25*up_samp_fac; % 25; 
dt = 16.5975 ./ 1000;
binsize = dt/up_samp_fac;
NT = NFRAMES*up_samp_fac;
[Ui_old, XVi_old] = NIM.generate_XVfolds( NT ); % original split of original data

% sanity check
subplot(1,3,1)
histogram(spk_bootstrap)
subplot(1,3,2)
histogram(spk_rnd) % somehow spk_rnd is better than spk_bootstrap?? maybe should not round(rate)
subplot(1,3,3)
spikes = datarun.spikes{3};
Robs = NIM.Spks2Robs(spikes, binsize, NT );
histogram(Robs)
set(gcf, 'Position', get(0, 'Screensize'));

%% set up training set vs test set of Xstim & fake spk
% [~,tt]=NIM.generate_XVfolds(20) % test set is always the third fifth of original data (at its middle)
mov = mov(:, :, XVi_old); % new mov is corresponding to fake data, aka old data's test set
[NY, NX, NFRAMES] = size(mov); NT = NFRAMES*up_samp_fac; % redefine using fake data size
[Ui, XVi] = NIM.generate_XVfolds( NT ); % re-split training vs test set after redefining NT
% [train_id, test_id] = NIM.generate_XVfolds( size(spikes, 1) );

RF = datarun.vision.sta_fits{cell_index,1};
RF_center_x = RF.mean(1);
RF_center_y = NY - RF.mean(2);
re_size = 15/2 - 1; % 15: taken from xml filename

x_upper = ceil(RF_center_x + re_size);
x_lower = floor(RF_center_x - re_size);
y_upper = ceil(RF_center_y + re_size);
y_lower = floor(RF_center_y - re_size);

mov_resize = mov(y_lower:y_upper, x_lower:x_upper, :);
NX = size(mov_resize,1); NY = NX;
mov = reshape(mov_resize, [NY*NX, NFRAMES])';
mov(mov < 0.5) = -0.48; % taken from xml contrast value
mov(mov > 0.5) = 0.48;

% Init stim parameters structure
params_stim = NIM.create_stim_params([nLags NY, NX], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);

% generate stimulus matrix
Xstim = NIM.create_time_embedding(mov,params_stim); % might run out of memory

%% round 1 fit
disp('Fitting GLM: 1 excitatory filter')
mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)
NL_types = {'lin'}; % define subunit as linear 

spikes = spk_bootstrap(1); % w less noise bc of bootstrapping, if there is no bias it should perform nicely
Robs = NIM.Spks2Robs(spikes, binsize, NT );

% New object-based definition of model (note regularization can be set here)
fit0_re = NIM( params_stim, NL_types, mod_signs, 'd2xt', 100 ); 
% also using spatial regularization too -- can do seprately or as one

% Tweak optimization parameters to doesn't creep around the minimum forever (necessary for 
% some reason with this data)
optim_params.maxIter = 400;

% Fit stimulus filters
fit0_re = fit0_re.fit_filters( Robs, Xstim, Ui, 'optim_params', optim_params );

% find optimal regularization using nested cross-val (note this can take a while)
% fit_test = fit0.reg_path( Robs, Xstim, Ui, XVi, 'lambdaID', 'd2xt', 'optim_params', optim_params );
% Subunit 1: Best reg = 600.00 for reduction in space mov; 400 in time; 200
% for spatial RF
% this finds optimal value of XX

% look at filter
k0 = reshape(fit0_re.subunits(1).filtK, [nLags, NY * NX]);
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

spikes = spk_bootstrap(2); 
Robs = NIM.Spks2Robs(spikes, binsize, NT);

% Try spike history term
fitS_re = fit0_re.init_spkhist(16,'doubling_time',4, 'negcon'); % I don't believe in positive spk-history
fitS_re = fitS_re.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );

spikes = spk_bootstrap(3);
Robs = NIM.Spks2Robs(spikes, binsize, NT);

% Always best to seed inhibition with different possibilities: random will often get stuck in local minima
fit1_re = fit0_re;
fit1_re.subunits(1).NLtype='rectlin'; % note new name (from threshlin)
fit1_re.subunits(2) = fit1_re.subunits(1); % make second subunit like the first
fit1_re.subunits(2).weight = -1; % make suppressive
kshift = NIM.shift_mat_zpad(k0,4,1); % shift filter to be delayed
fit1_re.subunits(2).filtK = kshift(:); % flatten

fit1_re = fit1_re.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
%fit1 = fit1.reg_path( Robs, Xstim, Ui, XVi, 'lambdaID', 'd2xt', 'fit_offsets', 1, 'optim_params', optim_params );

% Fit with spike-history term

spikes = spk_bootstrap(4); % expect fit1S still have highest R2 & LL?
Robs = NIM.Spks2Robs(spikes, binsize, NT);

fit1S_re = fit1_re.init_spkhist(16,'doubling_time',4, 'negcon');
fit1S_re = fit1S_re.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );
% spike history term helps
fit1S_re.display_model('Xstims', Xstim )
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit1S-' num2str(cell_id) '.png'])
close
% plot by itself
% plot(fit1S.spk_hist.bin_edges(2:end), fit1S.spk_hist.coefs)
% close

% new model display
fit1_re.display_model('Xstims',Xstim)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit1-' num2str(cell_id) '.png'])
close

% can also plot RFs separately
k1 = reshape(fit1_re.subunits(1).filtK, [nLags, NY*NX]);
k2 = reshape(fit1_re.subunits(2).filtK, [nLags, NY*NX]);
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

spikes = spk_bootstrap(1); % fit0nl = fit0? no spk_hist.coefs in model param
Robs = NIM.Spks2Robs(spikes, binsize, NT);
fit0nl_re = fit0_re.fit_spkNL( Robs, Xstim, Ui );

spikes = spk_bootstrap(3); 
Robs = NIM.Spks2Robs(spikes, binsize, NT);
fit1nl_re = fit1_re.fit_spkNL( Robs, Xstim, Ui );

% fit upstream nonlinearity?
spikes = spk_bootstrap(5); 
Robs = NIM.Spks2Robs(spikes, binsize, NT);

fit2_re = fit1_re.init_nonpar_NLs( Xstim );
fit2_re = fit2_re.fit_upstreamNLs( Robs, Xstim, Ui );
% this actually ends up looking rectlin, so not worth it
fit2_re.display_model('Xstims',Xstim)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit2-' num2str(cell_id) '.png'])
close

% internal LLx (cross-validated likelihood, although was used for meta-params)
[LLs,~,~,fp] = fit0_re.eval_model(Robs, Xstim, XVi );
LLs(2) = fitS_re.eval_model(Robs, Xstim, XVi );
LLs(3) = fit1_re.eval_model(Robs, Xstim, XVi );
LLs(4) = fit1S_re.eval_model(Robs, Xstim, XVi );
LLs(5) = fit2_re.eval_model(Robs, Xstim, XVi ); % fit2 (w nonlin) is not better than fit1
LLs-fp.nullLL

save after_refit_bootstrap_63.mat
