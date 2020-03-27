%% load data 
clear
clc

dataset_num = '/data000-map-sorted'; 
date_num = '2020-02-29-0';

% prefix_now = '/Volumes/dusom_fieldlab';
prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, dataset_num, dataset_num);
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);

%% load ds cell 

load('ds_cell_map_20200306.mat', 'ds_map_all');
flag = find(ds_map_all(:,1)==0);
ds_slave_id_mapPCA = ds_map_all(1:(flag(1)-1), 2); ismember(ds_slave_id_mapPCA, datarun.cell_ids)
ds_slave_id_mapEI = ds_map_all((flag(1)+1):(flag(2)-1), 2);
ds_slave_id_map2 = ds_map_all((flag(2)+1):end, 2); ds_slave_id_map2(ds_slave_id_map2 == 0) = [];
% ds_map_all

%% chop data000 into sections

gaps = round(diff(datarun.triggers));
switch_flag = gaps ~= 2 & gaps ~= 4;
tmp_triggers = datarun.triggers(switch_flag);
switch_index = [];
for i = 1 : length(tmp_triggers)
    switch_index = [switch_index, find(datarun.triggers == tmp_triggers(i))];
end
section_end = [300; datarun.triggers(switch_index); datarun.duration];
section_start = [0; datarun.triggers(1); datarun.triggers(switch_index + 1)];
sections = [section_start, section_end, (section_end - section_start)];
sections(end, 2:3) = [datarun.triggers(end)+4, datarun.triggers(end)+4 - sections(end, 1)];
sections(end+1, :) = [datarun.triggers(end)+4, datarun.duration, datarun.duration - (datarun.triggers(end)+4)];

ndf = [99; 5;5;5;5;5;5;5; 4;4;4;4;4; 3;3;3; 2; 99];
flash_config = [0; 2.2;2.8;2.4;2.1;2.2;2.4;2.8; 2.2;2.4;2.8;2.2;2.1; 2.2;2.4;4.8; 4.2; 0];
nflash = [sections(1:15,3)./2; sections(16:17,3)./4; sections(end,3)./2];
sections = [sections, ndf, flash_config, nflash];

sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);

%% 
slave_ds_id_all = unique(ds_map_all(:,2)); 
slave_ds_id_all(slave_ds_id_all == 0) = [];
% slave_ds_id_all

ds_slave_index = find(datarun.cell_ids == slave_ds_id_all(1)); 
spike_time = datarun.spikes{ds_slave_index, 1};

%%
binsize = 0.020;
binnum = datarun.duration / binsize;
edges = linspace(0, datarun.duration, binnum);
[binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each bin
% histcounts(binned,[0,1,2,3,4,5,6,7]) % see bins w 0-7 spikes within them

trial_len = 2 / binsize; % every trial = 2s = 100 bins
section_idx = [round(section_sort(:,1)/binsize), round(section_sort(:,2)/binsize), section_sort];

%%
trial_num = floor((section_idx(:,2) - section_idx(:,1)) / trial_len);
section_idx(:,end) = (section_idx(:,6) * 10 + section_idx(:,7));

% index sections w same flash intensity. cut off 2-4s

range_null = section_idx(find(section_idx(:,end)==990), 1:2)
marker = unique(section_idx(:,end), 'stable');

% range_flash = [];
% for i = 2 : length(marker)
%     range_flash = [range_flash; section_idx(find(section_idx(:,end)==marker(i)), 1:2)];
% end


