% 2018-09-26-0/ data007-nyj-map/ white noise
% data005 (195-4000s) & data006 flow 
clear
clc
close

%% convert data to spikes.mat & mov.mat
% 
% datapath = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/2018-09-26-0/data007-nyj-map/data007-nyj-map';
% datarun = load_data(datapath);
% datarun = load_neurons(datarun);
% datarun = load_params(datarun);
% 
% % test = cellfun(@length, datarun.spikes)
% % median(test) % since demo spikes length ~ 10^5, we'll try fit datarun.spikes{3}
% % mean(test)
% spikes = datarun.spikes{3};
% 
% movie_path = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/ds/code/NIM/BW-15-1-0.48-11111-53x40-60.35_xoffset2.xml'
% mvi = load_movie(movie_path, datarun.triggers);
% [mov, ~,~, dur, refresh] = get_movie_LL(movie_path, datarun.triggers, 216000); 
% % why? Nframe should theoretically be: datarun.duration * 60.35 = 217260
% 
% mov = squeeze(mov(:,:,1,:)); % no need for color dimension 
% [NY, NX, NFRAMES] = size(mov);
% mov = reshape(mov, [NX*NY, NFRAMES])';
% mov(mov < 0.5) = -0.48; % take from xml contrast value
% mov(mov > 0.5) = 0.48;
% 
% save('spikes_20180926_007_n3.mat', 'spikes') 
% save('mov_20180926_007.mat', 'mov', '-v7.3') % force save >2GB .mat

%% reload converted data
cd D:/RRR/Grad/Rotation/GF_lab/lab_Mac/ds/code/NIM
load spikes_20180926_007_n3 
load mov_20180926_007 

%%
up_samp_fac = 1; % 1; A little nicer resolution at 2, but indeed runs longer
tent_basis_spacing = 1;
nLags = 25*up_samp_fac; % 25; 
dt = 16.5975 ./ 1000;
binsize = dt/up_samp_fac;

% Stimulus is a T x M matrix of 1-d binary white noise (M bars-wide, T time steps)
[NFRAMES, Nstixel] = size(mov);
NT = NFRAMES*up_samp_fac;
NX = 53; NY = 40; % insert mov stim param: x=width, y=height

% Init stim parameters structure
params_stim = NIM.create_stim_params([nLags NX, NY], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);
% generate stimulus matrix
Xstim = NIM.create_time_embedding(mov,params_stim); % out of memory!
% bin spike times 
Robs = NIM.Spks2Robs(spikes, binsize, NT );

% For generating nested-cross validation indices (this is default 5-fold)
[Ui, XVi] = NIM.generate_XVfolds( NT ); % XVi array size = NT/5


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
fit_test = fit0.reg_path( Robs, Xstim, Ui, XVi, 'lambdaID', 'd2xt', 'optim_params', optim_params );
% this finds optimal value of XX

% look at filter
k0 = reshape(fit0.subunits(1).filtK, [nLags, NX*NY]);
[mval,bestlag] = max(max(abs(k0)')); bestlag


figure(2); clf;
[U, S, V] = svd(k0);
subplot(1,2,1)
plot(U(:,1))
axis tight
axis square
subplot(1,2,2); colormap gray
imagesc(reshape(V(:,1), NX, NY))
axis square
