
%% Setup
% Load data
load LGN_FFdata

% Set parameters of fits
up_samp_fac = 8; % temporal up-sampling factor applied to stimulus 
tent_basis_spacing = 2; % represent stimulus filters using tent-bases with this spacing (in up-sampled time units)
nLags = 48; % number of time lags for estimating stimulus filters

% Create structure with parameters using static NIM function
params_stim = NIM.create_stim_params( [nLags 1 1], 'stim_dt', DTstim, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing );

% Create T x nLags 'design matrix' representing the relevant stimulus history at each time point
Xstim = NIM.create_time_embedding( FFstim, params_stim );

% Bin spikes at analysis resolution (use internal function)
Robs = NIM.Spks2Robs( FFspks, params_stim.dt, size(Xstim,1) );
% equivalent operation: histc( FFspks, (0:(size(Xstim,1)-1))*params_stim.dt );

%% Fit a single-filter LN model (without cross-validation)
NL_types = {'lin'}; % define subunit as linear (note requires cell array of strings)
subunit_signs = [1]; % determines whether input is exc or sup (mult by +1 in the linear case)

% Set initial regularization as second temporal derivative of filter
lambda_d2t = 1;

% Initialize NIM 'object' (use 'help NIM.NIM' for more details about the contructor 
LN0 = NIM( params_stim, NL_types, subunit_signs, 'd2t', lambda_d2t );

% Fit model filters
LN0 = LN0.fit_filters( Robs, Xstim );

% Plot filter -- should be less smooth than what we would like (our 'prior')
filtfig = figure; hold on
plot(LN0.subunits(1).filtK,'r')

% Increase smoothness regularization and refit
%   you could do this: fit0.mods(1).reg_params.lambda_d2T = 80;
%   but we use the function NIM.set_reg_params:
LN0 = LN0.set_reg_params( 'd2t', 40 );
LN0 = LN0.fit_filters( Robs, Xstim );

figure(filtfig); hold on
plot(LN0.subunits(1).filtK,'k')

%% Cross-validation
% It is best to set aside some data for cross-validation, resulting in 'train_inds' and 'test_inds'
% Here, break off 1/5 of stimulus (in middle) to not use for fit
NT = length(Robs);
test_inds = ceil(NT*2/5):ceil(NT*3/5);
train_inds = setdiff(1:NT,test_inds);

LN0 = LN0.fit_filters( Robs, Xstim, train_inds );
% compare trained and tested log-likelihoods
[LN0.eval_model( Robs,Xstim, train_inds ) LN0.eval_model( Robs,Xstim, test_inds )]
% about the same, which is a good sign for not overfitting

% Now can use nested cross-validation to determine optimal regularization
LN1 = LN0.reg_path( Robs, Xstim, train_inds, test_inds );

%% GLM: add spike-history term and refit
GLM0 = LN1.init_spkhist( 20, 'doubling_time', 5 ); 
GLM0 = GLM0.fit_filters( Robs, Xstim, train_inds, 'silent', 1 ); % silent suppresses optimization output

% Display model components
GLM0.display_model('Xstims',Xstim,'Robs',Robs)

% Compare likelihoods
[LN0.eval_model( Robs, Xstim, test_inds )  GLM0.eval_model( Robs, Xstim, test_inds )]
% Also available in fit-structures
[LN0.fit_props.LL GLM0.fit_props.LL]
	
%% NIM: linear plus suppressive; like Butts et al (2011)
% Add an inhibitory input (with delayed copy of GLM filter and fit 
delayed_filt = NIM.shift_mat_zpad( LN0.subunits(1).filtK, 4 );
NIM0 = GLM0.add_subunits( {'rectlin'}, -1, 'init_filts', {delayed_filt} );
NIM0 = NIM0.fit_filters( Robs, Xstim, train_inds );

% Allow threshold of suppressive term to vary
NIM1 = NIM0.fit_filters( Robs, Xstim, 'fit_offsets', 1 ); % doesnt make huge difference
% Search for optimal regularization
NIM1 = NIM1.reg_path( Robs, Xstim, train_inds, test_inds, 'lambdaID', 'd2t' );

% Compare subunit filters (note delay of inhibition)
NIM1.display_subunit_filters()

%% Convert upstream NLs to nonparametric (monotonic) functions and fit
nonpar_reg = 20; % set regularization value
NIM2 = NIM1.init_nonpar_NLs( Xstim, 'lambda_nld2', nonpar_reg );
NIM2 = NIM2.fit_upstreamNLs( Robs, Xstim, 'silent', 1 );

% Do another iteration of fitting filters and upstream NLs
NIM2 = NIM2.fit_filters( Robs, Xstim, 'silent', 1 );
NIM2 = NIM2.fit_upstreamNLs( Robs, Xstim, 'silent', 1 );

% Optimal regularization
NIM2 = NIM2.reg_path( Robs, Xstim, train_inds, test_inds, 'lambdaID', 'nld2' );

%% Fit spiking Nonlinearity parameters
% the default model is a log-exp nonlinearity with Poisson likelihood 
NIM3 = NIM1.fit_spkNL( Robs, Xstim, train_inds );  % using parametric-upstreamNL model

% convert to Bernoulli with logistic spiking nonlinearity
NIM4 = NIM3;
NIM4.noise_dist = 'bernoulli';  % currently cannot fit spiking nonlinear with bernoulli likelihood
NIM4.spkNL.type = 'logistic';
NIM4 = NIM4.fit_filters( Robs, Xstim, train_inds );  % rescale filters accordingly
NIM4 = NIM4.fit_spkNL( Robs, Xstim, train_inds ); 

% For fun use quadratic model
GQM0 = GLM0.add_subunits( {'quad'}, -1, 'init_filts', {delayed_filt} ); % exc unit does much worse
GQM0 = GQM0.fit_filters( Robs, Xstim, train_inds );

% Optimize all models 
LNf = LN0.fit_spkNL( Robs, Xstim, train_inds, 'silent', 1 );  
GLMf = GLM0.reg_path( Robs, Xstim, train_inds, test_inds, 'lambdaID', 'd2t' );
GLMf = GLMf.fit_spkNL( Robs, Xstim, train_inds, 'silent', 1 );  
GQMf = GQM0.reg_path( Robs, Xstim, train_inds, test_inds, 'lambdaID', 'd2t' );
GQMf = GQMf.fit_spkNL( Robs, Xstim, train_inds, 'silent', 1 );  
NIMf = NIM3;

% Final cross-validation on nested data (this was used to fit meta-params)
[LNf.eval_model( Robs, Xstim, test_inds ) GLMf.eval_model( Robs, Xstim, test_inds ) ...
	GQMf.eval_model( Robs, Xstim, test_inds ) NIMf.eval_model( Robs, Xstim, test_inds )]

% Cross-validate on repeat data
XstimReps = NIM.create_time_embedding( FFstimR, params_stim );
RobsReps = NIM.Spks2Robs_reps( FFspksR, params_stim.dt, size(XstimReps,1) );

[LLs_LNM,LLnulls] = LNf.eval_model_reps( RobsReps, XstimReps );
LLs_GLM = GLMf.eval_model_reps( RobsReps, XstimReps );
LLs_GQM = GQMf.eval_model_reps( RobsReps, XstimReps );
LLs_NIM = NIMf.eval_model_reps( RobsReps, XstimReps );
LLs_NIMb = NIM4.eval_model_reps( RobsReps, XstimReps );

% plot cross-validated likelihoods across trials relative to 'null' model
figure; hold on
plot(LLs_LNM-LLnulls,'r')
plot(LLs_GLM-LLnulls,'c')
plot(LLs_GQM-LLnulls,'b')
plot(LLs_NIM-LLnulls,'g','LineWidth',1)
plot(LLs_NIMb-LLnulls,'g--','LineWidth',1)
legend('LNM','GLM','GQM','NIM')
% Note that the first few trials are worse for NIM and GQM, as described in Butts et al (2011)


%% Finally, look at PSTH comparisons to calculate R-squared
% Simulate and compare predicted firing rates
dt = params_stim.dt;
T = length(FFstimR)*DTstim;
[~,~,ratesLNM] = LNf.eval_model_reps( RobsReps, XstimReps );
[~,~,ratesGLM] = GLMf.eval_model_reps( RobsReps, XstimReps );
[~,~,ratesNIM] = NIMf.eval_model_reps( RobsReps, XstimReps );
%[~,~,ratesNIMb] = NIM4.eval_model_reps( RobsReps, XstimReps );

Nreps = size(ratesLNM,2);
psthLNM = sum(ratesLNM,2)/Nreps/dt;
psthGLM = sum(ratesGLM,2)/Nreps/dt;
psthNIM = sum(ratesNIM,2)/Nreps/dt;
%psthNIMb = sum(ratesNIMb,2)/Nreps/dt;
psthOBS = sum(RobsReps,2)/Nreps/dt;

ts = 1070:1200;
figure; hold on
plot(ts*dt,psthOBS(ts),'k','LineWidth',1);
plot(ts*dt,psthLNM(ts),'r');
plot(ts*dt,psthGLM(ts),'c');
plot(ts*dt,psthNIM(ts),'g');
%plot(ts*dt,psthNIMb(ts),'g--');

%% R2
R2s(1) = 1-mean((psthOBS-psthLNM).^2)/var(psthOBS);
R2s(2) = 1-mean((psthOBS-psthGLM).^2)/var(psthOBS);
R2s(3) = 1-mean((psthOBS-psthNIM).^2)/var(psthOBS)
%R2s(4) = 1-mean((psthOBS-psthNIMb).^2)/var(psthOBS)
