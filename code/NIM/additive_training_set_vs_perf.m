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


% reduce Xstim size in time
re_time = 216000 / 60 * 10; % take first 10 mins
mov_resize = mov(:, :, 1:re_time);
mov = reshape(mov_resize, [NY*NX, re_time])';


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
