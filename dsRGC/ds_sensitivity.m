% sensitivity plot for single cell
% later will import ds cell id list to achieve sequential plotting for
% every single ds cell

%% load data w changeable dataset_num & ds_now
clear
clc

dataset_num = '00-map'; % dim flashes to test absolute sensitivity
prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data0', dataset_num, '/data0', dataset_num);
% datapath = '/Volumes/???/lab/Experiments/Array/Analysis/2019-11-21-0/data002/data002';
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);

%% load ds cell identified in master 

load('ds_cell_map_20200207.mat', 'ds_cells', 'ds_map_all');
flag = find(ds_map_all(:,1)==0);
ds_slave_id_mapPCA = ds_map_all(1:(flag(1)-1), 2);
ds_slave_id_mapEI = ds_map_all((flag(1)+1):(flag(2)-1), 2);
ds_slave_id_map2 = ds_map_all((flag(2)+1):end, 2);

%% select cell & chop data000 into sections

ds_slave_id = ds_slave_index_map2(1);
ds_slave_index = find(datarun_s.cell_ids == ds_slave_id);
spike_time = datarun.spikes{ds_slave_index, 1};

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
sections = [sections(1:12,:); [sections(12,2), sections(13,1), (sections(13,1)-sections(12,2))]; sections(13:end,:)];

ndf = [99; 5;4;5; 4;4;4;4;4; 3;4;3;99; 3;2;2;99];
flash_config = [0; 2.4;2.4;2.8; 2.2;2.4;2.2;2.8;2.4; 2.2;2.8;2.4;0; 2.8;4.2;4.8;0];
nflash = [0; sections(2:12,3)./2; 0; sections(14,3)./2; sections(15:16,3)./4; 0];
sections = [sections, ndf, flash_config, nflash];

sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);
marker = unique(section_sort(:,7));
marker_seq = section_sort(:,7);

%% merge sections w same NDF and flash_config. x_axis=2 after cutting off 2-4s

for m = 1: (length(marker))
    subplot(length(marker), 1, m)  
    marker_now = marker(m);
    section_id_seq = find(marker_seq == marker_now, length(marker_seq));
    
    for s = 1:length(section_id_seq)
        section_id = section_id_seq(s); 
        section_now = [section_sort(section_id, 1), section_sort(section_id, 2)];
        section_flag = spike_time >= section_now(1) & spike_time <= section_now(2);
        
        spike_time_section = spike_time(section_flag);
        spike_time_section = spike_time_section - section_now(1);

        rep_len = 2;
        rep_max = round(section_now(2) - section_now(1)) / rep_len;
        
        if section_sort(section_id, 5) < 3  % period = 2s or 0s (dark)
            rep_step = 1;
        elseif section_sort(section_id, 5) >= 4 % period = 4s
            rep_step = 2; % skip 2-4s of every period, plot only 0-2s
        end
        
        for rep = 1 : rep_step : rep_max
            rep_flag = spike_time_section >= (rep-1)*rep_len & spike_time_section <= rep*rep_len;
            spike_time_rep = spike_time_section(rep_flag);
            spike_time_rep = spike_time_rep - (rep-1)*rep_len;

            rep_mark = rep * ones(length(spike_time_rep),1);
            scatter(spike_time_rep, rep_mark, 10, 'filled')
            axis([-0.05 (rep_len + 0.05) 0 (rep_max + 1)])
            hold on
        end
        hold on
    end
end

