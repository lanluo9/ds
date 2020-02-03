%% load data w changeable dataset_num & ds_now
clear
clc

% dataset_num = '00'; % dim flashes to test absolute sensitivity
prefix_now = '/Volumes/dusom_fieldlab';

% datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data0', dataset_num, '/data0', dataset_num);
% datapath = '/Volumes/???/lab/Experiments/Array/Analysis/2019-11-21-0/data002/data002';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/000-map-002-200129/000-map-002-200129');
% /Volumes/dusom_fieldlab/lab/Experiments/Array/Analysis/2019-11-21-0/000-map-002-200129

datarun = load_data(datapath);
datarun = load_neurons(datarun);
% datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);

%%
% cell_master_id_mapped = [469 2867 3710 4399 4621 5105 6318 6423 6695 7291]; % result of map-analysis
cell_master_id_mapped = [469 2869 3710 4399 4621 5105 6320 6423 6695 7291]; % shifted bc master002 was spike sorted
% [sanity, check] = ismembertol(cell_master_id_mapped(2), datarun.cell_ids, 1)
[sanity, check] = ismembertol(2869, datarun.cell_ids, 1)