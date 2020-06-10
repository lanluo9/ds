%% Make sure the NIMtoolbox is in your matlab path
% Load data: stimulus, spike times, and time res of stimulus
clear all
load RGCsimdata.mat

%% Format data and model parameters

% Format data for fitting
nLags = 50; % number of time lags for estimating stimulus filters
up_samp_fac = 8; % temporal up-sampling factor applied to stimulus 
tent_basis_spacing = 3; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
params_stim = NIMcreate_stim_params( nLags, stim_dt, up_samp_fac, tent_basis_spacing );

% Create T x nLags 'Xmatrix' representing the relevant stimulus history at each time point
Xstim = create_time_embedding( stim, params_stim ); 

% Bin spikes at analysis resolution
Robs = histc(spiketimes,(0:(size(Xstim,1)-1))*params_stim.dt); % Prob w conventions

%% Create and fit a (regularized) GLM (single linear filter LN model)
mod_signs = [1]; % determines whether input is exc or sup (doesn't matter in the linear case)
NL_types = {'lin'}; % define subunit as linear 

% Create regularization parameter structure (in this case we're only using a
% smoothness penalty in time (2nd deriv with respect to the time)
params_reg = NIMcreate_reg_params( 'lambda_d2T', 200 ); 
silent = 0; % set to 1 if you want to suppress the optimization display

% Initialize NIM (use 'help NIMinitialize_model' for more details about the 
% model struct components).
fit0 = NIMinitialize_model( params_stim, mod_signs, NL_types, params_reg ); 

% Fit stimulus filters
fit0 = NIMfit_filters( fit0, Robs, Xstim, [], [], silent ); 

% Display linear filter (first one) -- SHOULD WE BE USING "MODS"
filtfig = figure;
plot(fit0.mods(1).filtK,'k');

%% Fit an NIM with two rectified excitatory inputs
mod_signs = [1 1]; % both inputs are excitatory. (+1 for excitatory, -1 for suppressive)
NL_types = {'threshlin','threshlin'}; % make both filters have threshold-linear upstream NLs

% Initialize model and fit stimulus filters -- random initialization
fit1 = NIMinitialize_model( params_stim, mod_signs, NL_types, params_reg );
fit1 = NIMfit_filters( fit1, Robs, Xstim, [], [], silent );

figure(filtfig); hold on; 
plot(fit1.mods(1).filtK,'b')
plot(fit1.mods(2).filtK,'r')

% Built-in function for displaying model components
NIMdisplay_model( fit1 )
% Add 'Xstim' to see how nonlinearities compare with stim distributions
NIMdisplay_model( fit1, Xstim )


%% Convert upstream NLs to nonparametric (monotonic) functions and fit
to_nonpar = [1 2]; % which subunits to convert to tent-basis (nonparametric) upstream NLs
lambda_NLd2 = 20; % smoothness regularization on the tent-basis coefs
fit2 = NIMinitialize_upstreamNLs( fit1, Xstim, to_nonpar, lambda_NLd2 ); % initializes the tent-basis representation

% Now fit the NLs (+ spk history term) -- where is spike history term?
fit2 = NIMfit_upstreamNLs( fit2, Robs, Xstim, [],[],[], silent ); % fit the upstream NL

%% Do another iteration of fitting filters and upstream NLs
silent = 1; % switching to stealth mode
fit3 = NIMfit_filters( fit2, Robs, Xstim, [],[], silent );

fit3 = NIMfit_upstreamNLs( fit3, Robs, Xstim, [],[],[], silent );

% Note that LLs aren't improving much with these additional optimizations
fprintf('Log-likelihoods so far: \n');
disp(fit3.LL_seq);

%% Fit spiking NL params
% Note that in this case re-estimating the spiking NL shape doesnt improve things 
% although this is not always the case.
spk_NL_display = 1; % turn on display to get a before-and-after picture of the spk NL
fit4 = NIMfit_logexp_spkNL( fit3, Robs, Xstim, [], spk_NL_display );

%% Fit quadratic model (with linear and squared terms)
mod_signs = [1 1]; % both excitatory (although doesnt matter for linear)
NL_types = {'lin','quad'}; 

% Initialize model and fit stimulus filters
quad0 = NIMinitialize_model( params_stim, mod_signs, NL_types, params_reg );
quad0 = NIMfit_filters( quad0, Robs, Xstim, [],[], silent );

%% View different models with NIMdisplay_model

NIMdisplay_model( fit4, Xstim ); %Plot NIM
NIMdisplay_model( quad0, Xstim ); %Plot quad model

%% Compare model fits directly with the true filters
est_filt_tax = (0:tent_basis_spacing:(nLags-1)*tent_basis_spacing)*params_stim.dt;
true_filt_tax = (0:(size(true_filts,2)-1))*sim_dt;

% These are the true filters
figure
subplot(2,1,1);hold on
plot(true_filt_tax,true_filts(1,:));
plot(true_filt_tax,true_filts(2,:),'r');
xlim([0 1])
line([0 1],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',12)
title('True filters');

% These are the NIM filters
subplot(2,1,2);hold on
Kmat = [fit4.mods(:).filtK];
minvals = min(Kmat);
[~,off_filt_id] = min(minvals);
on_filt_id = setdiff(1:2,off_filt_id);
plot(est_filt_tax,fit4.mods(on_filt_id).filtK,'o-')
plot(est_filt_tax,fit4.mods(off_filt_id).filtK,'ro-')
xlim([0 1]);
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('Time lag (s)','fontsize',12)
title('NIM filter estimates')

% Repeat comparison for the quadratic model
figure
subplot(2,1,1);hold on
plot(true_filt_tax,true_filts(1,:));
plot(true_filt_tax,true_filts(2,:),'r');
xlim([0 1])
line([0 1],[0 0],'color','k','linestyle','--');
xlabel('Time lag (s)','fontsize',12)
title('True filters');
subplot(2,1,2);hold on
plot(est_filt_tax,quad0.mods(1).filtK,'o-')
plot(est_filt_tax,quad0.mods(2).filtK,'ro-')
xlim([0 1]);
line([0 1],[0 0],'color','k','linestyle','--')
xlabel('Time lag (s)','fontsize',12)
title('Quadratic model filter estimates')

%%