datapath = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-11-21-0/data006/data006';
%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun.names.stimulus_path = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s06.txt';

%% process triggers and extract some stim params
trigger_set = round(datarun.triggers);
trig_inds = find(mod(trigger_set, 10) == 0);
user_defined_trigs = datarun.triggers(trig_inds);
datarun = load_stim(datarun, 'user_defined_trigger_set', trig_inds);
num_stim = length(datarun.stimulus.combinations);
grat_dirs = datarun.stimulus.params.DIRECTION;
grat_TPs = datarun.stimulus.params.TEMPORAL_PERIOD;
stim_dur = 8;
num_reps = datarun.stimulus.repetitions;
num_rgcs = length(datarun.cell_ids);

%%
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

scatter((vector_mags_120), (vector_mags_240))


%%
% set x-y cuoff
cutoff_coord = [1.2, 1.2]; % leads to 47 failed_to_map out of 49 dsRGC in data006
cutoff_coord = [1.0, 1.0]; % leads to 50 failed_to_map out of 53 dsRGC in data006

cutoff_coord = [0.9, 1.2]; % leads to 25 failed_to_map out of 28 dsRGC in data002
cutoff_coord = [1.0, 1.0]; % leads to 29 failed_to_map out of 32 dsRGC in data002


x_finder = find(vector_mags_120 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun.cell_ids(selected_indices);


%%
slave_path = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-11-21-0/data000/data000';
datarun_s = load_data(slave_path);
datarun_s = load_neurons(datarun_s);
datarun_s = load_params(datarun_s);
datarun_s = load_ei(datarun_s, 'all', 'array_type', 519);

%%

[map_list, failed_to_map_list] = map_ei(datarun, datarun_s, 'master_cell_type', ds_cell_ids, 'slave_cell_type', 'all');
length(failed_to_map_list)
length(ds_cell_ids)




