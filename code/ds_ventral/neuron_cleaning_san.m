%% changeable dataset_num
clear
clc

dataset_num = '02-sorted';
date_num = '2020-02-29-0';

prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/dusom_fieldlab/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/data0', dataset_num, ...
    '/data0', dataset_num); 

stimulus_path = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/stimuli/s', ...
    dataset_num, '.txt');

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
% datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);

%% process triggers and extract some stim params
trigger_set = round(datarun.triggers);
trig_inds = find(mod(trigger_set, 10) == 0);
user_defined_trigs = datarun.triggers(trig_inds);

datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_set', trig_inds);

num_stim = length(datarun.stimulus.combinations);
grat_dirs = datarun.stimulus.params.DIRECTION;
grat_TPs = datarun.stimulus.params.TEMPORAL_PERIOD;
stim_dur = 8;
num_reps = datarun.stimulus.repetitions;
num_rgcs = length(datarun.cell_ids);

g_sp = datarun.stimulus.params.SPATIAL_PERIOD(1);
for g_dirs = 1:length(datarun.stimulus.params.DIRECTION)
    for g_tp = 1:length(datarun.stimulus.params.TEMPORAL_PERIOD)
        tmp_dir = datarun.stimulus.params.DIRECTION(g_dirs);
        tmp_tp = datarun.stimulus.params.TEMPORAL_PERIOD(g_tp);
        temp_spike_times = get_grating_spike_times(datarun, 'all', stim_dur,...
            'direction', tmp_dir, 'TP', tmp_tp, 'SP', g_sp);
        gratingrun.direction(g_dirs).temporal_period(g_tp).spike_times = temp_spike_times;
    end
end

%% scatter plot of vector sums for two different gratings.
[vector_sums_120, vector_mags_120] = get_vector_sums(datarun, 'all', 'TP', 120, 'SP', 240);
[vector_sums_240, vector_mags_240] = get_vector_sums(datarun, 'all', 'TP', 240, 'SP', 240);

% scatter((vector_mags_120), (vector_mags_240))

%%
x_cutoff = 0.7;
y_cutoff = 0.8;
cutoff_coord = [x_cutoff, y_cutoff]; 

close
scatter((vector_mags_120), (vector_mags_240))
hold on
xline(x_cutoff);
yline(y_cutoff);

x_finder = find(vector_mags_120 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun.cell_ids(selected_indices);
ds_index = selected_indices;
ds_cells = [ds_index; ds_cell_ids];
ds_cells'

%%
ds_master_id_mapPCA = [212; 2027; 2612; 3108; 3257; 3318; 3588; 3664; 3707; 3948; 4381; 5133; 5898; 6841; 7442];

for i = 1 : length(ds_master_id_mapPCA)
    single_ds_id = ds_master_id_mapPCA(i); 
    single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
    
    if ~isempty(single_ds_index)
        disp([num2str(single_ds_id), '  found in datarun.cell_id'])
    else
        disp([num2str(single_ds_id), '  NOT found in datarun.cell_id!!!'])        
    end
end
