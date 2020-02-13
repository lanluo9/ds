clear
clc
% contamination < 0.1
% spike rate > 0.01

strprovide = struct('load_neurons',1,'load_params',0,'load_ei',0);

datarun = load_data('/Volumes/dusom_fieldlab/lab/Experiments/Array/Analysis/2019-11-21-0/rerun_san/data002-sorted/data002-sorted',strprovide); 
datarun_s = load_data('/Volumes/dusom_fieldlab/lab/Experiments/Array/Analysis/2019-11-21-0/rerun_san/data000-map/data000-map',strprovide); 
id_kept = intersect(datarun.cell_ids, datarun_s.cell_ids);

ds_cells_002_sorted = load('ds_master_002_sorted_20200211.mat'); 
% load('ds_cell_map_20200210.mat', 'ds_map_all'); 

ds_cell_ids = ds_cells_002_sorted.ds_cells(2,:);
ds_kept = intersect(id_kept, ds_cell_ids);
ds_lost = setdiff(ds_cell_ids,intersect(id_kept, ds_cell_ids));

length(intersect(id_kept, ds_cell_ids))
length(ds_cell_ids)

% length(datarun.cell_ids)
% length(datarun_s.cell_ids)
% length(id_kept)


%%
clear
% contamination < 100
% spike rate > 0.00
strprovide = struct('load_neurons',1,'load_params',0,'load_ei',0);

datarun = load_data('/Volumes/dusom_fieldlab/lab/Experiments/Array/Analysis/2019-11-21-0/rerun_san/data002-sorted/data002-sorted',strprovide); 
datarun_ss = load_data('/Volumes/dusom_fieldlab/lab/Experiments/Array/Analysis/2019-11-21-0/rerun/data000-map/data000-map',strprovide); 
id_kept = intersect(datarun.cell_ids, datarun_ss.cell_ids);

ds_cells_002_sorted = load('ds_master_002_sorted_20200211.mat'); 
ds_cell_ids = ds_cells_002_sorted.ds_cells(2,:);
intersect(id_kept, ds_cell_ids)
length(intersect(id_kept, ds_cell_ids))
length(ds_cell_ids)

% length(datarun.cell_ids)
% length(datarun_ss.cell_ids)
% length(id_kept)
