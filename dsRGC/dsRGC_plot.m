%% changeable dataset_num & ds_now
clear
clc

dataset_num = '02';
prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data0', dataset_num, '/data0', dataset_num);
% datapath = '/Volumes/???/lab/Experiments/Array/Analysis/2019-11-21-0/data002/data002';

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun.names.stimulus_path = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s', dataset_num, '.txt');
% datarun.names.stimulus_path = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s02.txt';

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
% cutoff_coord = [1.2, 1.2]; % leads to 47 failed_to_map out of 49 dsRGC in data006
% cutoff_coord = [1.0, 1.0]; % leads to 50 failed_to_map out of 53 dsRGC in data006

% cutoff_coord = [0.9, 1.2]; % leads to 25 failed_to_map out of 28 dsRGC in data002
cutoff_coord = [1.0, 1.0]; % leads to 29 failed_to_map out of 32 dsRGC in data002

x_finder = find(vector_mags_120 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun.cell_ids(selected_indices);

%%
slave_path = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data000/data000');
datarun_s = load_data(slave_path);
datarun_s = load_neurons(datarun_s);
datarun_s = load_params(datarun_s);
datarun_s = load_ei(datarun_s, 'all', 'array_type', 519);

%%
[map_list, failed_to_map_list] = map_ei(datarun, datarun_s, 'master_cell_type', ds_cell_ids, 'slave_cell_type', 'all', 'troubleshoot', true);
fprintf('failed to map %d neurons out of %d neurons \n', length(failed_to_map_list), length(ds_cell_ids)); 

% %% dsRGC tuning curve Cartesian w merged TPs
% 
% dir_spike_count = [];
% for dir = 1:length(grat_dirs)
%     for tp = 1:length(grat_TPs)
%         tp_spike_count = 0;
%         single_dirtp = gratingrun.direction(dir).temporal_period(tp).spike_times;
%         neuron_spike_count = zeros(size(single_dirtp,1),1);
%         for neuron = 1:size(single_dirtp,1)
%             for rep = 1:size(single_dirtp,2)
%                 neuron_spike_count(neuron,1) = neuron_spike_count(neuron,1) + ...
%                     length(gratingrun.direction(dir).temporal_period(tp).spike_times{neuron,rep});
%             end
%         end
%         tp_spike_count = tp_spike_count + neuron_spike_count;
%     end
%     dir_spike_count = [dir_spike_count, tp_spike_count];
% end
% 
% ds_index = selected_indices;
% ds_spike_count = [];
% for i = 1:length(ds_index)
%     ds_spike_count = [ds_spike_count; dir_spike_count(ds_index(i),:)];
% end
% 
% for i = 1:length(ds_index)
%     plot(grat_dirs, ds_spike_count(i,:))
%     hold on
% end
% 
% %% dsRGC tuning curve polar plot w merged TPs
% 
% theta = deg2rad(grat_dirs);
% theta = [theta, theta(1)];
% % radius = ds_spike_count(1,:);
% % radius = [radius, radius(1)];
% % polarplot(theta, radius)
% 
% for i = 1:length(ds_index)
%     radius = ds_spike_count(i,:);
%     radius = [radius, radius(1)];
%     polarplot(theta, radius)
%     hold on
% end

%% 

ds_now = 30; % range 1:length(ds_index)
tp_set = 2; % range 1:length(grat_TPs), in this case 1:3

ds_index = selected_indices;
ds_cells = [ds_index; ds_cell_ids];
single_ds_index = ds_cells(1,ds_now);
single_ds_id = ds_cells(2,ds_now);

%%
savefile = 'ds_cells.mat';
save(savefile, 'ds_cells');

%% rasterplot by direction for single dsRGC w separate TPs

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

ds_spike_count = [];
for i = 1:length(ds_index)
    ds_spike_count = [ds_spike_count; dir_spike_count(ds_index(i),:)];
end

subplot(3,3,5)
theta = deg2rad(grat_dirs);
theta = [theta, theta(1)];
radius = ds_spike_count(find(ds_index==single_ds_index),:);
radius = [radius, radius(1)];
polarplot(theta, radius)
title(['data0', num2str(dataset_num), '. dsRGC index = ', num2str(single_ds_index), '. id = ', num2str(single_ds_id), '. TP = ', num2str(tp_set)])

