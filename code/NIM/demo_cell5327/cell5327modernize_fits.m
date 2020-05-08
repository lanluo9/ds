%% Modernized Field-Retina fit example

% cd ~/Desktop/cell5327/
load spikes
load mov

up_samp_fac = 1; % 1; A little nicer resolution at 2, but indeed runs longer
tent_basis_spacing = 1;
nLags = 25*up_samp_fac; % 25; 
dt = 16.5975 ./ 1000;
binsize = dt/up_samp_fac;

% Stimulus is a T x M matrix of 1-d binary white noise (M bars-wide, T time steps)
[NFRAMES, nXPix] = size(mov);
NT = NFRAMES*up_samp_fac;
NX = 11;  % can handle 2-d stim, just need to specify in stim params:

% Init stim parameters structure
params_stim = NIM.create_stim_params([nLags NX, NX], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);
%params_stim = NIMcreate_stim_params([nLags nXPix],dt, up_samp_fac, tent_basis_spacing);

% generate stimulus matrix
Xstim = NIM.create_time_embedding(mov,params_stim); 
%Xstim = create_time_embedding(mov,params_stim); 

% bin spike times 
spikes = spike{3,1};
Robs = NIM.Spks2Robs(spikes, binsize, NT );

% Could shift several time steps (6) forward, to get rid of early zeros in RF
% (but will leave as-is for now, for clarity. NIM.shift_mat_zpad is an easy function...

% For generating nested-cross validation indices (this is default 5-fold)
[Ui, XVi] = NIM.generate_XVfolds( NT );


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
%fit_test = fit0.reg_path( Robs, Xstim, Ui, XVi, 'lambdaID', 'd2xt', 'optim_params', optim_params );
% this finds optimal value of XX

% look at filter
k0 = reshape(fit0.subunits(1).filtK, [nLags, NX*NX]);
[mval,bestlag] = max(max(abs(k0)')); bestlag


figure(2); clf;
[U, S, V] = svd(k0);
subplot(1,2,1)
plot(U(:,1))
axis tight
axis square
subplot(1,2,2); colormap gray
imagesc(reshape(V(:,1),NX,NX))
axis square

%% round 2 fit NIM
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
% plot by itself
figure
plot(fit1S.spk_hist.bin_edges(2:end), fit1S.spk_hist.coefs)




% new model display
fit1.display_model('Xstims',Xstim)
% can also plot RFs separately
k1 = reshape(fit1.subunits(1).filtK, [nLags, NX*NX]);
k2 = reshape(fit1.subunits(2).filtK, [nLags, NX*NX]);
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
imagesc(reshape(k1(blat1,:),NX,NX),[-1 1]*mk1)
axis square
subplot(1,3,3); colormap gray
imagesc(reshape(k2(blat2,:),NX,NX), [-1 1]*mk2)
axis square

fit0nl = fit0.fit_spkNL( Robs, Xstim, Ui );
fit1nl = fit1.fit_spkNL( Robs, Xstim, Ui );

% fit upstream nonlinearity?
fit2 = fit1.init_nonpar_NLs( Xstim );
fit2 = fit2.fit_upstreamNLs( Robs, Xstim, Ui );
% this actually ends up looking rectlin, so not worth it
fit2.display_model('Xstims',Xstim)

% save models
save Res2mod3.mat fit0 fit1

% internal LLx (cross-validated likelihood, although was used for meta-params)
[LLs,~,~,fp] = fit0.eval_model(Robs, Xstim, XVi );
LLs(2) = fitS.eval_model(Robs, Xstim, XVi );
LLs(3) = fit1.eval_model(Robs, Xstim, XVi );
LLs(4) = fit1S.eval_model(Robs, Xstim, XVi );
LLs-fp.nullLL



%% Cross-validation
load temp_raster
load mov_rpts

Nreps = length(temp_raster);
% make raster, my style
spks_obs = [];
for rr = 1:Nreps
	spks_obs = [spks_obs; temp_raster{rr}; -1;];
end
NTR = size(mov_rpts,1);
% Bin to make Robs -- note shift by latenct change
RobsR = NIM.shift_mat_zpad(	histc(spks_obs,(0:(NTR-1))*binsize)/Nreps, 4 );
% weird thing is that stimulus/data has to be shifted by 32 ms. Could there
% be an alignment problem with the repeats? (or did I do something dumb?)

% Generate model predictions
Xstim_rep = NIM.create_time_embedding(mov_rpts,params_stim); 
[~,Rpred0] = fit0.eval_model([], Xstim_rep );
[~,Rpred1] = fit1.eval_model([], Xstim_rep );


% Spike history terms must be evaluated using input spike trains
% (note I dont see a spike-history "simulator" (like in the first NIM, which I could mix up if useful)
% right now just makes predictions based on previous observed spikes....
Robs_rep = NIM.Spks2Robs_reps(spks_obs'+4*binsize, binsize, NTR );
[LLsR,LLsRnull,rep_rates] = fit1S.eval_model_reps( Robs_rep, Xstim_rep );
Rpred1S = mean(rep_rates')';

% align PSTHs
bg_rg_d = 50; % data selection
ed_rg_d = 300;
bg_rg_s = 48; % sim selection
ed_rg_s = 298;
R2s = 1-[var(RobsR(bg_rg:ed_rg)-Rpred0(bg_rg_s:ed_rg_s)) var(RobsR(bg_rg:ed_rg)-Rpred1(bg_rg_s:ed_rg_s)) var(RobsR(bg_rg:ed_rg)-Rpred1S(bg_rg_s:ed_rg_s))]/var(RobsR(bg_rg:ed_rg))


figure(5); clf; hold on
plot((1:NTR-50)*dt/2, RobsR(bg_rg:ed_rg)/binsize,'k','LineWidth',1)
plot((1:NTR-50)*dt/2, Rpred0(bg_rg_s:ed_rg_s)/binsize,'r')
plot((1:NTR-50)*dt/2, Rpred1(bg_rg_s:ed_rg_s)/binsize,'g')
plot((1:NTR-50)*dt/2, Rpred1S(bg_rg_s:ed_rg_s)/binsize,'b')

legend('data', 'GLM', 'NIM', 'NIM-spk')
xlim([1.3 2.6]) % particularly good section
