%% changeable dataset_num
clear
clc

dataset_num = '02-sorted';

prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/dusom_fieldlab/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data0', dataset_num, ...
    '/data0', dataset_num); 

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun.names.stimulus_path = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s', ...
    dataset_num, '.txt');

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

%% test case
% cell index 8: low firing neuron w abnormal burst of high firing in 1 repetition (both TP 120 & 240)
% cell index 1: high firing neuron w abnormal uniformly added spikes in 1 direction (both TP 120 & 240)
% cell index 4: high firing neuron w abnormal uniformly added spikes in all directions (both TP 120 & 240)

datarun_fake = datarun;
i = 1;
while i <= 600
    datarun_fake.spikes{8,1}(end+1) = 0.28 + i/1000000;
    datarun_fake.spikes{8,1}(end+1) = 2.98 + i/1000000;
    i = i+1;
end
datarun_fake.spikes{8,1} = sort(datarun_fake.spikes{8,1});
datarun = datarun_fake;

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

%% test consistency code
% [vector_sums_120_original, vector_mags_120_original] = get_vector_sums(datarun, 'all', 'TP', 120, 'SP', 240);
% get_vector_sums_similarity_index
[vector_sums_120, vector_mags_120] = get_vector_sums(datarun, 'all', 'TP', 120, 'SP', 240);
[vector_sums_240, vector_mags_240] = get_vector_sums(datarun, 'all', 'TP', 240, 'SP', 240);

scatter((vector_mags_120), (vector_mags_240))

%%
% set x-y cuoff
cutoff_coord = [0.9, 0.9]; 

x_finder = find(vector_mags_120 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun.cell_ids(selected_indices);
ds_index = selected_indices;
ds_cells = [ds_index; ds_cell_ids];

% savefile = append('ds_master_002_sorted_', datestr(now, 'yyyymmdd'), '.mat');
% save(savefile, 'ds_cells');

%% rasterplot by direction for single dsRGC w separate TPs
% single_ds_id = ds_master_id_mapPCA(8); 
single_ds_id = datarun_fake.cell_ids(8);
tp_set = 1; % range 1:length(grat_TPs), in this case 1:3
% single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
single_ds_index = 8;

dir_spike_count = [];
for dir = 1:length(grat_dirs)
    for tp = tp_set
        tp_spike_count = 0;
        single_dirtp = gratingrun.direction(dir).temporal_period(tp_set).spike_times;
        neuron_spike_count = zeros(size(single_dirtp,1),1);
        for neuron = 1:size(single_dirtp,1)
            for rep = 1:size(single_dirtp,2)
                neuron_spike_count(neuron,1) = neuron_spike_count(neuron,1) + ...
                    length(gratingrun.direction(dir).temporal_period(tp_set).spike_times{neuron,rep});
            end
        end
        tp_spike_count = tp_spike_count + neuron_spike_count;
    end
    dir_spike_count = [dir_spike_count, tp_spike_count];
end

%%

subplot_num = [6 3 2 1 4 7 8 9; grat_dirs];
for dir = 1:length(grat_dirs) 
    subplot(3,3,subplot_num(1,dir))
    single_dirtp = gratingrun.direction(dir).temporal_period(tp_set).spike_times;
    for rep = 1:size(single_dirtp,2)        
        spike_time = gratingrun.direction(dir).temporal_period(tp_set).spike_times{single_ds_index, rep};
        rep_mark = rep .* ones(length(spike_time),1);
        scatter(spike_time, rep_mark, 50, 'filled')
        axis([0 10 0 7])
        hold on
    end
end

subplot(3,3,5)
theta = deg2rad(grat_dirs);
theta = [theta, theta(1)];
% radius = ds_spike_count(find(ds_index==single_ds_index),:);
radius = dir_spike_count(single_ds_index,:);
radius = [radius, radius(1)];
polarplot(theta, radius)
title(['data0', num2str(dataset_num), '. dsRGC index = ', num2str(single_ds_index), '. id = ', num2str(single_ds_id), '. TP = ', num2str(tp_set)])

