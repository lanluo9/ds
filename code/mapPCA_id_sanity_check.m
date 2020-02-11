clear
clc

dataset_num = '02-sorted';
% dataset_num = '02';

prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';
% prefix_now = '/Volumes/dusom_fieldlab/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/rerun/data0', dataset_num, ...
    '/data0', dataset_num); % rerun
% datapath = '/Volumes/???/lab/Experiments/Array/Analysis/2019-11-21-0/data002/data002';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
% datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);
datarun.names.stimulus_path = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s', dataset_num, '.txt');
% datarun.names.stimulus_path = '/Volumes/All_Staff/lab/Experiments/Array/Analysis/2019-11-21-0/stimuli/s02.txt';

