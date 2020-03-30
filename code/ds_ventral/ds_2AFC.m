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
ds_slave_id_mapPCA = ds_map_all(1:(flag(1)-1), 2); % ismember(ds_slave_id_mapPCA, datarun.cell_ids)
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
section_end = [300; datarun.triggers(switch_index); (datarun.duration-30)];
section_start = [0; datarun.triggers(1); datarun.triggers(switch_index + 1)];
sections = [section_start, section_end, (section_end - section_start)];
sections(end, 2:3) = [datarun.triggers(end)+4, datarun.triggers(end)+4 - sections(end, 1)];
sections(end+1, :) = [datarun.triggers(end)+4, datarun.duration-30, datarun.duration-30 - (datarun.triggers(end)+4)];

ndf = [99; 5;5;5;5;5;5;5; 4;4;4;4;4; 3;3;3; 2; 99];
flash_config = [0; 2.2;2.8;2.4;2.1;2.2;2.4;2.8; 2.2;2.4;2.8;2.2;2.1; 2.2;2.4;4.8; 4.2; 0];
nflash = [sections(1:15,3)./2; sections(16:17,3)./4; sections(end,3)./2];
sections = [sections, ndf, flash_config, nflash];

sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);

%% for a single cell

slave_ds_id_all = unique(ds_map_all(:,2)); 
slave_ds_id_all(slave_ds_id_all == 0) = [];
% slave_ds_id_all
ds_slave_index = find(datarun.cell_ids == slave_ds_id_all(1)); 

spike_time = datarun.spikes{ds_slave_index, 1};
binsize = 0.020; % divide into 20 ms bins
binnum = datarun.duration / binsize;
edges = linspace(0, datarun.duration, binnum);
[binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each bin
% histcounts(binned,[0:1:max(binned)+1]) % see bins w certain number of spikes within them

section_idx = [round(section_sort(:,1)/binsize), round(section_sort(:,2)/binsize), section_sort];
section_idx(:,end) = (section_idx(:,6) * 10 + section_idx(:,7));

%% iterate flash intensities

% trial_len = unique(floor(section_idx(fid_seq(1),7))) / binsize;
trial_len = 2 / binsize; % for now, pretend all trial length = 2s

sum_null = zeros(ntrial(1), trial_len); 
for t = 1 : ntrial(1)
    trial_null = binned(trial_len*(t-1)+section_idx(1,1)+1 : trial_len*t+section_idx(1,1));
    sum_null(t,:) = trial_null;
end
sum_null_start = sum(sum_null,1);
sum_null = zeros(ntrial(2), trial_len); 
for t = 1 : ntrial(2)
    trial_null = binned(trial_len*(t-1)+section_idx(2,1)+1 : trial_len*t+section_idx(2,1));
    sum_null(t,:) = trial_null;
end
sum_null_end = sum(sum_null,1);

ntrial = round(section_idx(:,8));
marker = unique(section_idx(:,end), 'stable');
Pc = zeros(length(marker)-1 ,1);
for flash_intensity = 2 : length(marker) % exclude dark==990
    fid = section_idx(:,end)==marker(flash_intensity);
    fid_seq = find(fid==1);
    if section_idx(fid_seq(1),6) == 5
        nid = 1; sum_null = sum_null_start;
    else
        nid = 2; sum_null = sum_null_end;
    end
    
    sum_flash_seq = zeros(1, trial_len);
    for i = 1 : length(fid_seq)
        sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len);
        for t = 1 : ntrial(fid_seq(i))
            trial_flash = binned(trial_len*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*t+section_idx(fid_seq(i),1));
            sum_flash_section{i}(t,:) = trial_flash;
        end
        sum_flash_section{i} = sum(sum_flash_section{i},1);
        sum_flash_seq = sum_flash_seq + sum_flash_section{i};
    end
    
    trial_num_null = 1 : ntrial(nid);
    trial_num_flash = 1 : sum(ntrial(fid));
    sample_size = min(length(trial_num_null), length(trial_num_flash));
    order_null = datasample(trial_num_null, sample_size, 'Replace', false);
    order_flash = datasample(trial_num_flash, sample_size, 'Replace', false);
    
    corrpos = zeros(sample_size, 1);
    for t = 1 : sample_size
        trial_null = binned(trial_len*(order_null(t)-1)+section_idx(nid,1)+1 : trial_len*order_null(t)+section_idx(nid,1));
        other_null = sum_null - trial_null;
        mean_null = other_null ./ (length(trial_num_null) - 1);

        if order_flash(t) <= ntrial(fid_seq(1))
            trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                trial_len*order_flash(t)+section_idx(fid_seq(1),1));
        else
            trial_flash = binned(trial_len*(order_flash(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                trial_len*(order_flash(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
        end
        other_flash = sum_flash_seq - trial_flash;
        mean_flash = other_flash ./ (length(trial_num_flash) - 1);

        discriminant = (mean_flash - mean_null)';
        corrpos(t) = (trial_flash - trial_null) * discriminant;
    end
    corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
    Pc(flash_intensity - 1) = corr / length(corrpos);
    
end

Pc

