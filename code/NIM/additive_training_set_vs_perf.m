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

cell_index = 3; %cell_index_array(counter);
cell_id = datarun.cell_ids(cell_index)
spikes = datarun.spikes{cell_index};

%%
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
clear mov_resize

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
[Ui_old, XVi_old] = NIM.generate_XVfolds( NT ); % XVi array size = NT/5. XVi is test set, while Ui is training set

%%
frame_per_sec = NFRAMES / datarun.duration; 
test_set_len = 10 * 60 * frame_per_sec; % fix test set length to 10 min
training_set_len = [1, 5, 10, 20, 30, 40, 50]  * 60 * frame_per_sec;

XVi = [(NT - test_set_len + 1) : 1 : NT]';
% Ui_seq_rnd = {};
Ui_seq_last = {};
for i = 1 : length(training_set_len)
    fin_max = min(XVi) - 1;
%     bgn = randi([1, fin_max - training_set_len(i) + 1]);    
%     Ui_seq_rnd{i,1} = [bgn : 1 : (bgn + training_set_len(i) - 1)];
    Ui_seq_last{i,1} = [(fin_max - training_set_len(i) + 1) : 1 : fin_max];
%     intersect(Ui_seq_last{i,1}, XVi)
end

%%
for j = 1 : length(training_set_len)
Ui = Ui_seq_last{j,1};

training_min = training_set_len(j) / frame_per_sec / 60

% %% round 1 fit
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
saveas(gcf, ['k0_filter-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
close

% %% round 2 fit NIM
mod_signs = [1 -1]; % both inputs are excitatory. (+1 for excitatory, -1 for suppressive)
NL_types = {'rectlin','rectlin'}; % name-change

% Try spike history term
fit0S = fit0.init_spkhist(16,'doubling_time',4, 'negcon'); % I don't believe in positive spk-history
fit0S = fit0S.fit_filters( Robs, Xstim, Ui, 'fit_offsets', 1, 'optim_params', optim_params );

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
saveas(gcf, ['fit1S-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
close
% plot by itself
% plot(fit1S.spk_hist.bin_edges(2:end), fit1S.spk_hist.coefs)
% close


% new model display
fit1.display_model('Xstims',Xstim)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit1-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
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
saveas(gcf, ['subunits-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
close

fit0nl = fit0.fit_spkNL( Robs, Xstim, Ui );
fit1nl = fit1.fit_spkNL( Robs, Xstim, Ui );
fit0Snl = fit0S.fit_spkNL( Robs, Xstim, Ui );
fit1Snl = fit1S.fit_spkNL( Robs, Xstim, Ui );

% fit upstream nonlinearity?
fit2 = fit1.init_nonpar_NLs( Xstim );
fit2 = fit2.fit_upstreamNLs( Robs, Xstim, Ui );
% this actually ends up looking rectlin, so not worth it
fit2.display_model('Xstims',Xstim)
set(gcf, 'Position', get(0, 'Screensize'));
saveas(gcf, ['fit2-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
close

% internal LLx (cross-validated likelihood, although was used for meta-params)
[LLs,~,~,fp] = fit0.eval_model(Robs, Xstim, XVi );
LLs(2) = fit0S.eval_model(Robs, Xstim, XVi );
LLs(3) = fit0nl.eval_model(Robs, Xstim, XVi );
LLs(4) = fit0Snl.eval_model(Robs, Xstim, XVi );

LLs(5) = fit1.eval_model(Robs, Xstim, XVi );
LLs(6) = fit1S.eval_model(Robs, Xstim, XVi );
LLs(7) = fit1nl.eval_model(Robs, Xstim, XVi );
LLs(8) = fit1Snl.eval_model(Robs, Xstim, XVi );

LLs(9) = fit2.eval_model(Robs, Xstim, XVi ); % fit2 (w nonlin) is not better than fit1
LLfit = LLs - fp.nullLL

save(['fit_' num2str(training_min) 'min_cell63_RFsize13.mat'], '-regexp', 'fit*')
save(['fit_' num2str(training_min) 'min_cell63_RFsize13.mat'], 'XVi', 'LLfit', '-append')
disp(['finished training set len ' num2str(training_min) 'min' ])

end


%% visualize fit_spkNL

for j = 1 : length(training_set_len)

    training_min = training_set_len(j) / frame_per_sec / 60
    load(['fit_' num2str(training_min) 'min_cell63_RFsize13.mat'], 'fit0nl', 'fit1nl')
        
%     figure('units','normalized','outerposition',[0 0 1 1/2]) 
    fit0nl.display_model('Xstims', Xstim )
    set(gcf, 'Position', [0 0 15000 250]);
    saveas(gcf, ['fit0nl-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
    close
    
    fit1nl.display_model('Xstims', Xstim )
    set(gcf, 'Position', get(0, 'Screensize'));
    saveas(gcf, ['fit1nl-' num2str(training_min) 'min_' num2str(cell_id) '.png'])
    close

end
fprintf('done')

%% test LL error for 50 min

[LLs_test,~,~,fp_test] = fit0.eval_model(Robs, Xstim, XVi_old );
LLs_test(2) = fit0S.eval_model(Robs, Xstim, XVi_old );
LLs_test(3) = fit0nl.eval_model(Robs, Xstim, XVi_old );
LLs_test(4) = fit0Snl.eval_model(Robs, Xstim, XVi_old );

LLs_test(5) = fit1.eval_model(Robs, Xstim, XVi_old );
LLs_test(6) = fit1S.eval_model(Robs, Xstim, XVi_old );
LLs_test(7) = fit1nl.eval_model(Robs, Xstim, XVi_old );
LLs_test(8) = fit1Snl.eval_model(Robs, Xstim, XVi_old );

LLs_test(9) = fit2.eval_model(Robs, Xstim, XVi_old ); % fit2 (w nonlin) is not better than fit1
LLfit_test = LLs_test - fp_test.nullLL

%%
[fit1.eval_model( Robs,Xstim, Ui_seq_last{end,1} ) fit1.eval_model( Robs,Xstim, XVi)]
