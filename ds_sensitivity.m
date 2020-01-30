% sensitivity plot for single cell
% later will import ds cell id list to achieve sequential plotting for
% every single ds cell

%% changeable dataset_num & ds_now
clear
clc

dataset_num = '00'; % dim flashes stimulus to test absolute sensitivity
prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data0', dataset_num, '/data0', dataset_num);
% datapath = '/Volumes/???/lab/Experiments/Array/Analysis/2019-11-21-0/data002/data002';

%% load data
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);
% datarun.names.stimulus_path = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s', dataset_num, '.txt');
% datarun.names.stimulus_path = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s02.txt';

%% process triggers and extract some stim params
% trigger_set = round(datarun.triggers);
% trig_inds = find(mod(trigger_set, 10) == 0);
% user_defined_trigs = datarun.triggers(trig_inds);
% datarun = load_stim(datarun, 'user_defined_trigger_set', trig_inds);
% num_stim = length(datarun.stimulus.combinations);
% grat_dirs = datarun.stimulus.params.DIRECTION;
% grat_TPs = datarun.stimulus.params.TEMPORAL_PERIOD;
% stim_dur = 8;
% num_reps = datarun.stimulus.repetitions;
% num_rgcs = length(datarun.cell_ids);

% verify datarun.spikes
% histogram(datarun.spikes{1,1}, 40)

%% scaffold: single section of darkness. cell id = 469, index = 21
gaps = round(diff(datarun.triggers));
switch_flag = gaps ~= 2 & gaps ~= 4;
switch_duration = [round(datarun.triggers(1)) - 300; gaps(switch_flag)];

tmp_triggers = datarun.triggers(switch_flag);
switch_index = [];
for i = 1 : length(tmp_triggers)
    switch_index = [switch_index, find(datarun.triggers == tmp_triggers(i))];
end
section_end_prev = [300; datarun.triggers(switch_index)];
section_start_next = [datarun.triggers(1); datarun.triggers(switch_index + 1)];
sections = [section_end_prev, section_start_next];

tolerance = 0.01;

%% scaffold: sanity check
% ntrig_flash = 60;
% nflash = [60];
% gap_seq = [300];
% for i = 2 : length(switch_index)
%     ntrig_flash = ntrig_flash + switch_index(i) - switch_index(i-1) -1;
%     nflash = [nflash, (switch_index(i) - switch_index(i-1) - 1)];
%     gap_seq = [gap_seq, (datarun.triggers(switch_index(i)) - datarun.triggers(switch_index(i-1)))];
% end
% % ntrig_flash + length(switch_index) + 1 == length(datarun.triggers)
% gap_seq = [gap_seq, (datarun.triggers(end) - gap_seq(end))];
% section_seq = transpose(gap_seq) - switch_duration;

