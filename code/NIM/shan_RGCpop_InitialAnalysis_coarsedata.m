%% Field Lab multi-neuron model fits -- initial demonstration

% cd ~/Desktop/
load coarse_white_noise_data_NDF0_4types_2018_11_30.mat
load indices_for_cpl_2018_11_30

% Establish experimental parameters
NC = length(indicies_for_cpl{1,1}{1,1});%细胞的总个数
LLstore = zeros(NC,5);
for cc = 4:5
NC_detail = indicies_for_cpl{1,1}{1,1};

white_noise_spikes_couple = white_noise_spikes{1,1}(NC_detail);%pick the connect cells

white_noise_movie = reshape(white_noise_movie,10,13,[]);
[NX, NY, NT] = size(white_noise_movie);%刺激噪声的三个维度
% organize stim to be +/- 1
stim = reshape(white_noise_movie,[NX*NY,NT])';%刺激矩阵，把刺激噪声先转换为:时间*空间
stim(stim < 0.5) = -1;
stim(stim > 0.5) = 1;

% % Bin spike train
frac = 1;
tbins = (0:frac*NT)*refresh/frac+white_noise_movie_start;%把仿真点换算成时间并加仿真开始时间，记录了每个刺激的开始时间
Robs = zeros(frac*NT, NC);%初始化放电为243000*10的空矩阵
% % Binned at double resolution
% % tbins2 = (0:2*NT)*refresh/2+white_noise_movie_start;
% % Robs2 = zeros(2*NT, NC); % second version at double-resolution (8 ms)
% % For 8x resolution
% frac = 8;
% Robs8 = zeros(frac*NT, NC); % second version at double-resolution (2 ms)
% tbins8 = (0:frac*NT)*refresh/frac+white_noise_movie_start;

for nn = 1:NC%一共10个神经元，这里做了循环写入
	dum = histc(white_noise_spikes_couple{nn}, tbins);%第1-10个细胞的放电按照每个刺激开始时间进行分段求和
	Robs(:,nn) = dum(1:end-1);%所有细胞放电求和后放入一个矩阵中（10列，每列为一个细胞）
% 	dum = histc(white_noise_spikes{nn}, tbins2);
% 	Robs2(:,nn) = dum(1:end-1);
% 	dum = histc(white_noise_spikes{nn}, tbins8);
% 	Robs8(:,nn) = dum(1:end-1);
end

% Save processed data
%save data_proc.mat stim NX NY Robs Robs2 Robs8 refresh

%% Model setup
lat_skip = 2;  % how many time lags to throw out preceding stimulus预测刺激之前的时滞
num_lags = 24; % number of time lags in temporal kernel: note that response is analyzed here at refresh/2
stim_par = NIM.create_stim_params( [num_lags NX NY], 'stim_dt', refresh, 'upsampling', frac );

% Xstim = NIM.create_time_embedding( stim, stim_par ); % design stim matrix
% Rsh = NIM.shift_mat_zpad( Robs2, -lat_skip );  % shifted responses
% [Ui, Xi] = NIM.generate_XVfolds( NT*2 ); % establish 5-fold cross-validation indices


Xstim = NIM.create_time_embedding( stim, stim_par ); % design stim matrix 
Rsh = NIM.shift_mat_zpad( Robs, -lat_skip );  % shifted responses 把放电往下移动2个时间点（lat）
[Ui, Xi] = NIM.generate_XVfolds( NT* frac); % establish 5-fold cross-validation indices 划分训练数据和验证数据


% STAs
stas = Xstim'*Rsh;
% figure
% for ccplot = 1:NC
% 	subplot(3,4, ccplot)
% 	k = reshape(stas(:,ccplot), [num_lags, NX*NY]);
% 	plot(k)
% end

% GLM
% cc = 3;  % fit example cell

glm0 = NIM( stim_par, {'lin'}, 1, 'd2xt', 100 );
glm0 = glm0.fit_filters( Rsh(:,cc), Xstim, Ui );%Rsh(:,cc)为第cc个神经元的放电
glm0 = glm0.reg_path( Rsh(:,cc), Xstim, Ui, Xi, 'lambdaID', 'd2xt' );

% Add spike-history term
glm1 = glm0.init_spkhist( 12, 'doubling_time', 4, 'negcon');
glm1 = glm1.fit_filters( Rsh(:,cc), Xstim, Ui );

% NIM (just playing here)
nim0 = NIM( stim_par, repmat({'rectlin'},2,1), [1 -1], 'd2xt', 100 );%
k0 = nim0.subunits(1).filtK;
nim0.subunits(1).filtK = k0;
nim0.subunits(2).filtK = -k0;
ksh = NIM.shift_mat_zpad(reshape(k0, [num_lags, NX*NY]), 3);
nim0.subunits(2).filtK = reshape(ksh, [num_lags*NX*NY, 1]);
nim0 = nim0.fit_filters( Rsh(:,cc), Xstim, Ui );
% nim0 = nim0.reg_path( Rsh(:,cc), Xstim, Ui, Xi, 'lambdaID', 'd2xt' ); % to optimize regularization

nim0b = nim0.fit_filters( Rsh(:,cc), Xstim, Ui, 'fit_offsets', 1 );
% nim0b.display_model('Xstims', Xstim)
% Fitting offsets had good effect on main excitatory kernel

% Add spike-history term
nim1 = nim0b.init_spkhist( 12, 'doubling_time', 4, 'negcon');
nim1 = nim1.fit_filters( Rsh(:,cc), Xstim, Ui, 'fit_offsets', 1 );

%% Adding time-lagged spike trains from other cells (causal)
% just demonstrating GLM here
other_ccs = setdiff(1:NC, cc); % other neurons 其他细胞的所在的列

num_coupling_lags = 24;  % will be at lower resolution其他细胞的输出作为刺激  lag
coupling_par = NIM.create_stim_params( [num_coupling_lags NC-1 1], 'stim_dt', refresh, 'tent_spacing', 2 );%其他神经元输出
Cstim = NIM.create_time_embedding( Robs(:,other_ccs), coupling_par );%
% makes a time-lagged version of the other spike trains

% Assemble design matrices
Xs{1} = Xstim;
Xs{2} =  Cstim;

% make a coupled GLM model
Pglm1 = glm1; 
Pglm1.stim_params(2) = coupling_par;  % add second stim params
Pglm1 = Pglm1.add_subunits( {'lin'}, [1], 'xtargs', 2, 'd2t', 10 );
Pglm1 = Pglm1.fit_filters( Rsh(:,cc), Xs, Ui );
% Pglm1.display_model('Xstims', Xs)

% make a couple NIM model (with auto-spike feedback)
Pnim = nim1; 
Pnim.stim_params(2) = coupling_par;  % add second stim params
Pnim = Pnim.add_subunits( {'lin'}, [1], 'xtargs', 2, 'd2t', 10 );
Pnim = Pnim.fit_filters( Rsh(:,cc), Xs, Ui );
% Pnim.display_model('Xstims', Xs)

% The big thing here is that a second stim-params was added, which corresponds to a different 'Xtarget'
% Each subunit can be specified to act on a different stim based on its setting 'xtarg' (i.e. Pglm1.subunits(1).xtarg)
% the model can be designed and fit from scratch as well, i.e.
% >> Pglm1 = NIM( [stim_par coupling_par], {'lin', 'lin'}, [1 1], 'Xtargets', [1 2], ...)

%% Model performance
[LLs,~,~,LLinfo] = glm0.eval_model(Rsh(:,cc), Xstim, Xi ); LLnull = LLinfo.nullLL;
LLs(2) = glm1.eval_model(Rsh(:,cc), Xstim, Xi );
LLs(3) = nim1.eval_model(Rsh(:,cc), Xstim, Xi );
LLs(4) = Pglm1.eval_model(Rsh(:,cc), Xs, Xi );
LLs(5) = Pnim.eval_model(Rsh(:,cc), Xs, Xi );
LLs-LLnull
LLstore(cc,:) = LLs-LLnull;
eval(['save D:\Matlabtest\offbt\OFFBTcell' num2str(cc) '.mat']); %change save path here
clearvars -except NC LLstore indicies_for_cpl movie_wnr refresh starts white_noise_movie white_noise_movie_start white_noise_spikes wnr_spikes
end


for i = 1:13
  eval(['LL' '=[LL' num2str(i)  '];']);
  row = size(LL,1);
  GLMdiff(1:row,i) = LL(:,4) - LL(:,2);
  NIMdiff(1:row,i) = LL(:,5) - LL(:,3);
  meanGLM(i) = sum(GLMdiff(1:row,i))/row;
  meanNIM(i) = sum(NIMdiff(1:row,i))/row;
end
x = 1:1:13;
plot (x, meanGLM)