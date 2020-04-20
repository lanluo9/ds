%% changeable dataset_num
clear
clc
close

% dataset_num = '02-sorted'; date_num = '2020-02-29-0';
% prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff';
% prefix_now = '/dusom_fieldlab/All_Staff';
% datapath = [prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/data0', dataset_num, '/data0', dataset_num]
% stimulus_path = [prefix_now, '/lab/Experiments/Array/Analysis/', date_num, '/stimuli/s02-sorted.txt'];

% cd \\duhsnas-pri.dhe.duke.edu\dusom_fieldlab\All_Staff\lab\Experiments\Array\Analysis\2020-02-29-0\data002-sorted
datapath = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/2020-02-29-0/data002-sorted/data002-sorted';
stimulus_path = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/2020-02-29-0/stimuli/s02-sorted.txt';

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);

trigger_set = round(datarun.triggers);
trig_inds = find(mod(trigger_set, 10) == 0);
user_defined_trigs = datarun.triggers(trig_inds);
datarun.names.stimulus_path = stimulus_path;
datarun = load_stim_gdf(datarun, 'user_defined_trigger_set', trig_inds);

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

% scatter((vector_mags_120), (vector_mags_240))
% hold on
% xline(x_cutoff);
% yline(y_cutoff);

x_finder = find(vector_mags_120 > cutoff_coord(1));
y_finder = find(vector_mags_240 > cutoff_coord(2));
selected_indices = intersect(x_finder, y_finder);

ds_cell_ids = datarun.cell_ids(selected_indices);
ds_index = selected_indices;
ds_cells = [ds_index; ds_cell_ids];
ds_cells'

%% load latency txt
normal_latency = importdata('latency.txt')
divide = find(normal_latency(:,1)==0);
ds_master_normal = normal_latency(1:(divide-1), 1);
ds_master_latency = normal_latency((divide+1):end, 1);
ds_master_all = [ds_master_normal; ds_master_latency];

%% central sudoku plot superimposed

for i = 1 : length(ds_master_all)
    figure
    
    single_ds_id = ds_master_all(i); 
    single_ds_index = ds_cells(1, ds_cells(2,:)==single_ds_id);
    if isempty(single_ds_index)
        disp([num2str(single_ds_id), ' not found in datarun.cell_id'])
        close
        continue
    end
    if ismember(single_ds_id, ds_master_normal)
        title_flag = 'normal';
    elseif ismember(single_ds_id, ds_master_latency)
        title_flag = 'long latency';
    end
    
    for j = 1:3
        tp_set = j;
        if tp_set == 1
            color = [1 0 0];
        elseif tp_set == 2
            color = [0 1 0];
        else 
            color = [0 0 1];
        end
        
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

        theta = deg2rad(grat_dirs);
        theta = [theta, theta(1)];
        radius = dir_spike_count(single_ds_index,:);
        radius = [radius, radius(1)];
        polarplot(theta, radius,'Color', color)
        hold on        
    end
    
    line([0, vector_sums_120(single_ds_index,1)], [0, vector_sums_120(single_ds_index,2)],'Color', color)
    line([0, vector_sums_240(single_ds_index,1)], [0, vector_sums_240(single_ds_index,2)],'Color', color)
    line([0, vector_sums_480(single_ds_index,1)], [0, vector_sums_480(single_ds_index,2)],'Color', color)
    
    title([title_flag ' ' num2str(single_ds_id)])
%     saveas(gcf, ['ds_across_TPs-', num2str(single_ds_id), '.jpg'])
%     print(['ds_across_TPs-', num2str(single_ds_id)], '-dpdf', '-fillpage')
    disp(['plotted fig for ', num2str(single_ds_id)])
%     close
end

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
%     title(['data0', num2str(dataset_num), '. dsRGC index = ', num2str(single_ds_index), '. id = ', num2str(single_ds_id), '. TP = ', num2str(tp_set)])
    
%     saveas(gcf, ['ds_across_TPs-', num2str(single_ds_id), '.jpg'])
%     print(['ds_across_TPs-', num2str(single_ds_id)], '-dpdf', '-fillpage')
%     disp(['saved fig for ', num2str(single_ds_id)])
%     close
end
