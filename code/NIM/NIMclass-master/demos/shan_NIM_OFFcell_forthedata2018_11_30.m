%% Make sure the NIMtoolbox is in your matlab path
% clear all
load coarse_white_noise_data_2018_11_30.mat
profile on
profile clear
%% pick cell
ind = 1; %(OFFBS:white_noise_spikes{1}:1-29;OFFBT:white_noise_spikes{2}:1-35)
celltype = 1 ;% Only select OFFBS=1 or OFFBT=2
%% dt

dt = refresh;
%% spike
if celltype == 1
    spikes_total = white_noise_spikes{1}{ind};
end
if celltype == 2
    spikes_total = white_noise_spikes{2}{ind};
end
%% choose some time for simulation
% slen = 100000;%just for simulation
% Selected_time = slen*dt;% simulation time (s)
% spikes = spikes_total(spikes_total<Selected_time);% spikes in simulation time
spikes = spikes_total - white_noise_movie_start;%white noise spikes have a slight timing offset with the white noise movie

%% movie
mov1 = reshape(white_noise_movie,10,13,[]);
mov2 = permute(mov1,[3 1 2]);% in NIM code , time first position

% %use pillow toolbox calculate STA 
% nkt = 30;
% Stim0 = reshape(mov1,10*13,[])';
% Stim = bsxfun(@minus,Stim0,mean(Stim0));
% sps_coarse=histc(spikes,(0:242999)*dt);
% sta = simpleSTC(Stim,sps_coarse,nkt);
% figure;imagesc(sta)
% %----------------------------------------------------

%% NIM.Format data and model parameters
% Format data for fitting
nLags = 30; % number of time lags for estimating stimulus filters
% up_samp_fac = 1; % temporal up-sampling factor applied to stimulus 
% tent_basis_spacing = 2; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
params_stim = NIM.create_stim_params( [nLags 10 13], 'stim_dt',dt ); 
% Create T x nLags 'Xmatrix' representing the relevant stimulus history at each time point
Xstim = NIM.create_time_embedding( mov2, params_stim ); 
% Bin spikes at analysis resolution 
Robs = histc(spikes,(0:(size(Xstim,1)-1))*params_stim.dt);
%% Modify the stimulus, Because this is an off-type neuron, This step will swap all 0.98 with 0.02
Xstim(Xstim>0.5)=-1;
Xstim(Xstim>0.01)=0.980000019073486;
Xstim(Xstim==-1)=0.020000010728836;

%% Add Cross-validation 
% It is best to set aside some data for cross-validation, resulting in 'train_inds' and 'test_inds'
% Here, break off 1/5 of stimulus (in middle) to not use for fit
NT = length(Robs);
% change train data to 75%
test_inds = ceil(NT*2/4):ceil(NT*3/4);%middle 1/5 used for test ,form 2/5---3/5 
train_inds = setdiff(1:NT,test_inds);%others used for train

%% 1.LN model:Create and fit a (regularized) GLM (single linear filter LN model) without spike-history
mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)
NL_types = {'lin'}; % define subunit as linear 
% Initialize
LNM = NIM( params_stim, NL_types, mod_signs, 'd2xt',10, 'l1',50 );
% Fit stimulus filters
LNM = LNM.fit_filters(Robs, Xstim, train_inds); 
% compare trained and tested log-likelihoods
[LNM.eval_model( Robs,Xstim, train_inds ) LNM.eval_model( Robs,Xstim, test_inds )]
% about the same, which is a good sign for not overfitting

% Now can use nested cross-validation to determine optimal regularization
LNM_reg = LNM.reg_path( Robs, Xstim, train_inds, test_inds );
LNM_reg.display_model('Xstims',Xstim,'Robs',Robs)
profile report
%% 2.GLM0: add spike-history term in LN and refit
GLM0 = LNM_reg.init_spkhist( 20, 'doubling_time', 5 ); 
GLM0 = GLM0.fit_filters( Robs, Xstim, train_inds); 
% GLM0.display_model('Xstims',Xstim,'Robs',Robs)
% Compare likelihoods
[LNM.eval_model( Robs, Xstim, test_inds )  GLM0.eval_model( Robs, Xstim, test_inds )]
% Also available in fit-structures
[LNM.fit_props.LL GLM0.fit_props.LL]
GLM_reg = GLM0.reg_path( Robs, Xstim, train_inds, test_inds );
GLM_reg.display_model('Xstims',Xstim,'Robs',Robs)

%% 3.Fit an NIM with two rectified excitatory inputs
mod_signs = [1 -1 ]; % both inputs are excitatory. (+1 for excitatory, -1 for suppressive)
NL_types = {'rectlin','rectlin'}; % make both filters have threshold-linear upstream NLs

% Initialize model and fit stimulus filters
NIM1 = NIM( params_stim, NL_types, mod_signs, 'd2xt',10, 'l1',50 );
NIM1 = NIM1.fit_filters(Robs, Xstim, train_inds );
% compare trained and tested log-likelihoods
[NIM1.eval_model( Robs,Xstim, train_inds ) NIM1.eval_model( Robs,Xstim, test_inds )]
% about the same, which is a good sign for not overfitting

% Use nested cross-validation to determine optimal regularization
NIM_reg = NIM1.reg_path( Robs, Xstim, train_inds, test_inds );

%% 4.add a spk history term to NIM
NIMaddhis = NIM_reg.init_spkhist( 20, 'doubling_time', 5 ); 
NIMaddhis = NIMaddhis.fit_filters( Robs, Xstim, train_inds); % silent suppresses optimization output
% NIMaddhis.display_model('Xstims',Xstim,'Robs',Robs)
% Compare likelihoods
[NIM1.eval_model( Robs, Xstim, test_inds )  NIMaddhis.eval_model( Robs, Xstim, test_inds )]
% Also available in fit-structures
[NIM1.fit_props.LL NIMaddhis.fit_props.LL]

NIMaddhis_reg = NIMaddhis.reg_path( Robs, Xstim, train_inds, test_inds );
% NIMaddhis_reg.display_model('Xstims',Xstim,'Robs',Robs)
[NIMaddhis.eval_model( Robs, Xstim, test_inds )  NIMaddhis_reg.eval_model( Robs, Xstim, test_inds )]

%% 5.Convert upstream NLs to nonparametric (monotonic) functions and fit
nonpar_reg = 20; % smoothness regularization on the tent-basis coefs
NIM3 = NIMaddhis_reg.init_nonpar_NLs( Xstim, 'lambda_nld2', nonpar_reg ); % initializes the tent-basis representation
NIM3 = NIM3.fit_upstreamNLs( Robs, Xstim); 
% Check result 
% NIM3.display_model('Xstims',Xstim,'Robs',Robs)

%% 6.Do another iteration of fitting filters and upstream NLs
NIM4 = NIM3.fit_filters( Robs, Xstim);
NIM4 = NIM4.fit_upstreamNLs( Robs, Xstim );
% Optimal regularization
NIM4 = NIM4.reg_path( Robs, Xstim, train_inds, test_inds, 'lambdaID', 'nld2' );
% Check result use NIM toolbox 
NIM4.display_model('Xstims',Xstim,'Robs',Robs)

%% 7.Fit spiking Nonlinearity parameters
% the default model is a log-exp nonlinearity with Poisson likelihood 
NIM5 = NIMaddhis_reg.fit_spkNL( Robs, Xstim, train_inds );  % using parametric-upstreamNL model
LNf = LNM.fit_spkNL( Robs, Xstim, train_inds, 'silent', 1 ); 
GLMf = GLM_reg.fit_spkNL( Robs, Xstim, train_inds, 'silent', 1 ); 
NIMf= NIM5;
% Final cross-validation on nested data (this was used to fit meta-params)
[LNf.eval_model( Robs, Xstim, test_inds ) GLMf.eval_model( Robs, Xstim, test_inds ) ...
	NIMf.eval_model( Robs, Xstim, test_inds )]


%% Cross-validate on repeat data
% get a cell's psth and compare it to model predictions
movRepeat = reshape(movie_wnr,300,10,13);
XstimReps = NIM.create_time_embedding( movRepeat, params_stim );
%% Because this is an off-type neuron, swap all 0.98 with 0.02 
XstimReps(XstimReps>0.5)=-1;
XstimReps(XstimReps>0.01)=0.980000019073486;
XstimReps(XstimReps==-1)=0.020000010728836;
%% find the spiking of repeated stimuli
%% spike_Repeat
if celltype == 1
    spikes_R_o = wnr_spikes{1}{ind};
end
if celltype == 2
    spikes_R_o = wnr_spikes{2}{ind};
end
spikes_R = spikes_R_o';
%% Processing data into NIN toolbox analyzable 
%(1)Divide the total length according to start and calculate from start time
one_trial_tltal_time = size(movie_wnr,1)* dt;

for n =1:length(starts)-1
    for m = 1:length(spikes_R)
        if spikes_R(m)> starts(n) && spikes_R(m) < starts(n)+one_trial_tltal_time
            spikes_R(m)=spikes_R(m)-starts(n);
        end
        if spikes_R(m)> starts(end)
            spikes_R(m)=spikes_R(m)-starts(end);
        end
    end
end
% (2)Insert -1 between every two experiments
for p = length(spikes_R):-1:2
   if spikes_R(p)<spikes_R(p-1)
       spikes_R(p+1:end+1)=spikes_R(p:end);
       spikes_R(p)=-1;
   end
end
spikes_R = [spikes_R -1];
%--------------------------Finished  Processing------------------------------
%% PSTH
RobsReps = NIM.Spks2Robs_reps( spikes_R, params_stim.dt, size(XstimReps,1) );
[~,~,ratesLNM] = LNf.eval_model_reps( RobsReps, XstimReps );
[~,~,ratesGLM] = GLMf.eval_model_reps( RobsReps, XstimReps );
[~,~,ratesNIM] = NIMf.eval_model_reps( RobsReps, XstimReps );
Nreps = size(ratesLNM,2);
psthLNM = sum(ratesLNM,2)/Nreps/dt;
psthGLM = sum(ratesGLM,2)/Nreps/dt;
psthNIM = sum(ratesNIM,2)/Nreps/dt;
psthOBS = sum(RobsReps,2)/Nreps/dt;
ts = 1:size(XstimReps,1);
figure; hold on
plot(ts*dt,psthOBS(ts),'k','LineWidth',1);
plot(ts*dt,psthLNM(ts),'c');
plot(ts*dt,psthGLM(ts),'r');
plot(ts*dt,psthNIM(ts),'g');
legend('OBS','LNM','GLM','NIM')
set(gca,'FontSize',14)
title(['OFFBT', ' cell', num2str(ind)])

%% Calculate five metrics compare models
PSTHOBS_2 = psthOBS(31:300)';%Cut off the 30 points at the beginning
PSTHLNM_2 = psthLNM(31:300)';
PSTHGLM_2 = psthGLM(31:300)';
PSTHNIM_2 = psthNIM(31:300)';
perfLNM = get_pred_fr_metrics( PSTHOBS_2,PSTHLNM_2);
perfGLM = get_pred_fr_metrics( PSTHOBS_2,PSTHGLM_2);
perfNIM = get_pred_fr_metrics( PSTHOBS_2,PSTHNIM_2);
likelihood = [LNf.eval_model( Robs, Xstim, test_inds ) GLMf.eval_model( Robs, Xstim, test_inds ) ...
	NIMf.eval_model( Robs, Xstim, test_inds )];
%¡ý£¨useless£©¡ýEasy data copying to document forms
record = zeros(5,3);
record(1:4,1)=perfLNM';
record(1:4,2)=perfGLM';
record(1:4,3)=perfNIM';
record(5,:)=likelihood ;
profile viewer