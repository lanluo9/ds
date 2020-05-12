% 2018-09-26-0/ data007-nyj-map/ white noise
% data005 (195-4000s) & data006 flow 
clear
clc
close

%% convert data to spikes.mat & mov.mat

datapath = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/2018-09-26-0/data007-nyj-map/data007-nyj-map';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);

% test = cellfun(@length, datarun.spikes)
% median(test) % since demo spikes length ~ 10^5, we'll try fit datarun.spikes{3}
% mean(test)
spikes = datarun.spikes{3};

movie_path = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/ds/code/NIM/BW-15-1-0.48-11111-53x40-60.35.xml'
mvi = load_movie(movie_path, datarun.triggers);
[mov, ~,~, dur, refresh] = get_movie_LL(movie_path, datarun.triggers, 216000); 
% why? Nframe should theoretically be: datarun.duration * 60.35 = 217260

mov = squeeze(mov(:,:,1,:)); % no need for color dimension 
[NX, NY, NFRAMES] = size(mov);
mov = reshape(mov, [NX*NY, NFRAMES])';
mov(mov < 0.5) = -0.48; % take from xml contrast value
mov(mov > 0.5) = 0.48;

%%
spikes_20180926_007_n3 = spikes;
save('spike.mat', 'spikes_20180926_007_n3') 
mov_20180926_007 = mov;
% save('nim_prep.mat', 'mov_20180926_007', '-append')
save('mov.mat', 'mov_20180926_007', '-v7.3')

%%
% find mov stim param
load spikes_2018_09_26_007 
load mov_2018_09_26_007 

up_samp_fac = 1; % 1; A little nicer resolution at 2, but indeed runs longer
tent_basis_spacing = 1;
nLags = 25*up_samp_fac; % 25; 
dt = 16.5975 ./ 1000;
binsize = dt/up_samp_fac;

% Stimulus is a T x M matrix of 1-d binary white noise (M bars-wide, T time steps)
[NFRAMES, nXPix] = size(mov);
NT = NFRAMES*up_samp_fac;
NX = 11;  % 2-d stim % insert mov stim param

% Init stim parameters structure
params_stim = NIM.create_stim_params([nLags NX, NX], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);
% generate stimulus matrix
Xstim = NIM.create_time_embedding(mov,params_stim); 
% bin spike times 
Robs = NIM.Spks2Robs(spikes, binsize, NT );
% add NIM.shift_mat_zpad to shift several time steps (6) forward, to get rid of early zeros in RF

% For generating nested-cross validation indices (this is default 5-fold)
[Ui, XVi] = NIM.generate_XVfolds( NT ); % XVi array size = NT/5

%%

