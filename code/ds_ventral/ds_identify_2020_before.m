%% changeable dataset_num
clear
clc

dataset_num = '02-sorted';
date_num = '2020-02-29-0';
% prefix_now = '/Volumes/dusom_fieldlab';
prefix_now = '/Volumes/All_Staff/';
% prefix_now = '/dusom_fieldlab/All_Staff/';

datapath = [prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/data0', dataset_num, '/data0', dataset_num]; 
stimulus_path = [prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/stimuli/s02-sorted.txt'];

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);

%% process triggers and extract some stim params
trigger_set = round(datarun.triggers);
trig_inds = find(mod(trigger_set, 10) == 0);
user_defined_trigs = datarun.triggers(trig_inds);

datarun.names.stimulus_path = stimulus_path;
datarun = load_stim(datarun, 'user_defined_trigger_set', trig_inds);

%%
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

[vector_sums_120, vector_mags_120] = get_vector_sums(datarun, 'all', 'TP', 120, 'SP', 240);
[vector_sums_240, vector_mags_240] = get_vector_sums(datarun, 'all', 'TP', 240, 'SP', 240);
[vector_sums_480, vector_mags_480] = get_vector_sums(datarun, 'all', 'TP', 480, 'SP', 240);

% scatter((vector_mags_120), (vector_mags_240))

%% scatter plot of vector sums for two different gratings.

x_cutoff = 0.8;
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
ds_cells = [ds_index; ds_cell_ids]';
% ds_cells
% savefile = append('ds_master_002_ds_', datestr(now, 'yyyymmdd'), '.mat');
% save(savefile, 'ds_cells');

%% load long latency vs normal cells (criteria: slave sensi rasterplot)
tmp = importdata('latency.txt');
ds_master_normal = tmp(1 : (find(tmp(:,1)==0)-1), 1);
ds_master_latency = tmp((find(tmp(:,1)==0)+1) : end, 1);

%% use vector_sums_120 & vector_sums_240 as preferred direction

pref_dir_120 = vector_sums_120(selected_indices,:);
pref_dir_240 = vector_sums_240(selected_indices,:);
pref_dir_480 = vector_sums_480(selected_indices,:);

figureHandle = figure;
subplot(1,3,1)
for i = 1 : length(selected_indices)
    line([0, pref_dir_120(i,1)], [0, pref_dir_120(i,2)])
    hold on
    xlim([-2,2])
    ylim([-2,2])
end
subplot(1,3,2)
for i = 1 : length(selected_indices)
    line([0, pref_dir_240(i,1)], [0, pref_dir_240(i,2)])
    hold on
    xlim([-1.5,2])
    ylim([-1.5,2])
end
subplot(1,3,3)
for i = 1 : length(selected_indices)
    line([0, pref_dir_480(i,1)], [0, pref_dir_480(i,2)])
    hold on
    xlim([-2,2])
    ylim([-2,2])
end
axesHandles = findobj(get(figureHandle,'Children'), 'flat','Type','axes'); % Extract axes handles of all subplots from figure
axis(axesHandles,'square') % Set axis property to square

%% speed tuning?

tp_axis = [120, 240, 480]; 
for i = 1 : length(selected_indices)
    plot(tp_axis, [vector_mags_120(selected_indices(i)), vector_mags_240(selected_indices(i)), vector_mags_480(selected_indices(i))])
    hold on
end

%%
slave_path = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/data000-map-sorted/data000-map-sorted');

datarun_s = load_data(slave_path);
datarun_s = load_neurons(datarun_s);
datarun_s = load_params(datarun_s);
datarun_s = load_ei(datarun_s, 'all', 'array_type', 519);

[map_list, failed_to_map_list] = map_ei_custom2(datarun, datarun_s, 'master_cell_type', ds_cell_ids, 'slave_cell_type', 'all', 'troubleshoot', true);
fprintf('failed to map %d neurons out of %d neurons \n', length(failed_to_map_list), length(ds_cell_ids)); 

%%
t = map_list(1:2,:)';
t2 = t(~cellfun('isempty', t));
ds_map_ei = cell2mat(reshape(t2,[length(t2)/2,2]));
ds_master_id_mapEI = ds_map_ei(:,1);                                      % 15 mapped by map_ei only, corr threshold 0.85
ds_slave_id_mapEI = ds_map_ei(:,2); 

% ds_master_id_mapPCA = importdata('mapPCA_id.txt');                      % 12 mapped by map-analysis only
ds_master_id_mapPCA = [2027; 2612; 3108; 3257; 3318; 3588; 3664; 3707; 3948; 4381; 6841; 7442]; 
ds_slave_id_mapPCA = ds_master_id_mapPCA; 
ismember(ds_master_id_mapPCA, datarun_s.cell_ids)

ds_master_id_map2 = intersect(ds_master_id_mapEI, ds_master_id_mapPCA);   % 10 survivors of both map-analysis & map_ei
ds_slave_id_map2 = intersect(ds_slave_id_mapEI, ds_slave_id_mapPCA);

patch = ismember(ds_master_id_map2, ds_slave_id_map2) .* ds_master_id_map2;
master_col = [ds_master_id_mapPCA; 0; ds_master_id_mapEI; 0; ds_master_id_map2];
slave_col = [ds_slave_id_mapPCA; 0; ds_slave_id_mapEI; 0; patch];
ds_map_all = [master_col, slave_col]

% savefile = append('ds_cell_map_', datestr(now, 'yyyymmdd'), '.mat');
% save(savefile, 'ds_cells', 'ds_map_all');

%% ds-ness of mapped cells


for i = 1 : length(ds_cell_ids)
    figure
    
    single_ds_id = ds_cell_ids(i); 
    single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
    if isempty(single_ds_index)
        disp([num2str(single_ds_id), ' not found in datarun.cell_id'])
        continue
    end
    tp_set = 2; 

    % should optimize this w cellfun to broadcast length function
    dir_spike_count = zeros(length(datarun.cell_ids), length(grat_dirs));
    for dir = 1:length(grat_dirs)
        single_dirtp = gratingrun.direction(dir).temporal_period(tp_set).spike_times;
        neuron_spike_count = zeros(size(single_dirtp,1),1);
        for neuron = 1:size(single_dirtp,1)
            for rep = 1:size(single_dirtp,2)
                neuron_spike_count(neuron,1) = neuron_spike_count(neuron,1) + ...
                    length(gratingrun.direction(dir).temporal_period(tp_set).spike_times{neuron,rep});
            end
        end
        dir_spike_count(:,dir) = neuron_spike_count;
    end

    subplot_num = [6 3 2 1 4 7 8 9; grat_dirs];
    for dir = 1:length(grat_dirs) 
        subplot(3,3,subplot_num(1,dir))
        for rep = 1:datarun.stimulus.repetitions        
            spike_time = gratingrun.direction(dir).temporal_period(tp_set).spike_times{single_ds_index, rep};
            rep_mark = rep .* ones(length(spike_time), 1);
            scatter(spike_time, rep_mark, 5, 'filled')
            axis([0 10 0 (datarun.stimulus.repetitions+1)])
            hold on
        end
    end

    subplot(3,3,5)
    theta = deg2rad(grat_dirs);
    theta = [theta, theta(1)];
    radius = dir_spike_count(single_ds_index,:);
    radius = [radius, radius(1)];
    polarplot(theta, radius)
    title(['data0', num2str(dataset_num), '. dsRGC index = ', num2str(single_ds_index), '. id = ', num2str(single_ds_id), '. TP = ', num2str(tp_set)])
    
    saveas(gcf, ['0229-', num2str(single_ds_index),'-', num2str(single_ds_id), '.png'])
    print(['0229-', num2str(single_ds_id)], '-dpdf', '-fillpage')
    disp(['saved fig for ', num2str(single_ds_id)])
    close
end

%%

for i = 1 : length(ds_master_id_mapEI)
    figure
    
    single_ds_id = ds_master_id_mapEI(i); 
    single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
    if isempty(single_ds_index)
        disp([num2str(single_ds_id), ' not found in datarun.cell_id'])
        continue
    end
    tp_set = 2; % range 1:3

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

    subplot_num = [6 3 2 1 4 7 8 9; grat_dirs];
    for dir = 1:length(grat_dirs) 
        subplot(3,3,subplot_num(1,dir))
        single_dirtp = gratingrun.direction(dir).temporal_period(tp_set).spike_times;
        for rep = 1:size(single_dirtp,2)        
            spike_time = gratingrun.direction(dir).temporal_period(tp_set).spike_times{single_ds_index, rep};
            rep_mark = rep .* ones(length(spike_time),1);
            scatter(spike_time, rep_mark, 5, 'filled')
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
    
%     saveas(gcf, ['mapEI-', num2str(single_ds_index),'-', num2str(single_ds_id), '.png'])
%     close
end

%%
for i = 1 : length(ds_master_id_mapPCA)
    figure
    
    single_ds_id = ds_master_id_mapPCA(i); 
    single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
    if isempty(single_ds_index)
        disp([num2str(single_ds_id), ' not found in datarun.cell_id'])
        continue
    end
    tp_set = 2; % range 1:3

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

    subplot_num = [6 3 2 1 4 7 8 9; grat_dirs];
    for dir = 1:length(grat_dirs) 
        subplot(3,3,subplot_num(1,dir))
        single_dirtp = gratingrun.direction(dir).temporal_period(tp_set).spike_times;
        for rep = 1:size(single_dirtp,2)        
            spike_time = gratingrun.direction(dir).temporal_period(tp_set).spike_times{single_ds_index, rep};
            rep_mark = rep .* ones(length(spike_time),1);
            scatter(spike_time, rep_mark, 5, 'filled')
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
    
%     saveas(gcf, ['mapPCA-', num2str(single_ds_index),'-', num2str(single_ds_id), '.png'])
%     close
end

%%
for i = 1 : length(ds_master_id_map2)
    figure
    
    single_ds_id = ds_master_id_map2(i); 
    single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
    if isempty(single_ds_index)
        disp([num2str(single_ds_id), ' not found in datarun.cell_id'])
        continue
    end
    tp_set = 2; % range 1:3

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

    subplot_num = [6 3 2 1 4 7 8 9; grat_dirs];
    for dir = 1:length(grat_dirs) 
        subplot(3,3,subplot_num(1,dir))
        single_dirtp = gratingrun.direction(dir).temporal_period(tp_set).spike_times;
        for rep = 1:size(single_dirtp,2)        
            spike_time = gratingrun.direction(dir).temporal_period(tp_set).spike_times{single_ds_index, rep};
            rep_mark = rep .* ones(length(spike_time),1);
            scatter(spike_time, rep_mark, 5, 'filled')
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
    
%     saveas(gcf, ['mapBoth-', num2str(single_ds_index),'-', num2str(single_ds_id), '.png'])
%     close
end