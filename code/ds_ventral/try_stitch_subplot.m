%% load data w changeable dataset_num & ds_now
clear
clc

dataset_num = '/data000-map'; 
date_num = '2016-03-04-0';

prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, dataset_num, dataset_num);
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);

%% chop data000 into sections

gaps = round(diff(datarun.triggers));
flash_period = 3;
switch_flag = gaps ~= flash_period;
tmp_triggers = datarun.triggers(switch_flag);
switch_index = [];
for i = 1 : length(tmp_triggers)
    switch_index = [switch_index, find(datarun.triggers == tmp_triggers(i))];
end
section_end = [300; datarun.triggers(switch_index); datarun.duration];
section_start = [0; datarun.triggers(1); datarun.triggers(switch_index + 1)];
sections = [section_start, section_end, (section_end - section_start)];
sections(end, 2:3) = [datarun.triggers(end)+flash_period, datarun.triggers(end)+flash_period - sections(end, 1)];
% sections(end+1, :) = [datarun.triggers(end)+flash_period, datarun.duration, datarun.duration - (datarun.triggers(end)+flash_period)];

ndf = [99; 5;5;5; 4;4;4;4;4; 3;3;3; 2;2; 1;0];
flash_config = [0; 3.2;3.4;3.8; 3.2;3.3;3.4;3.6;3.8; 3.2;3.4;3.8; 3.2;3.4; 3.2;3.2];
nflash = [0; sections(2:end,3)./flash_period];
sections = [sections, ndf, flash_config, nflash];

sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);
% marker = unique(section_sort(:,7));
% marker_seq = section_sort(:,7);

load('ds_master_002_ds_20200311.mat')
slave_ds_id_all = unique(ds_cells(2,:)); 

%%
tic

for i = 1 : length(slave_ds_id_all)
    figure('units','normalized','outerposition',[0 0 1 1]) 

    ds_slave_id = slave_ds_id_all(i); 
    ds_slave_index = find(datarun.cell_ids == ds_slave_id); 
    if isempty(ds_slave_index)
        disp([num2str(ds_slave_id), ' not found in slave datarun.cell_id'])
        close
        continue
    end
    
    spike_time = datarun.spikes{ds_slave_index, 1};

    for section_id = 1:size(section_sort, 1)
%         subplot(size(section_sort, 1), 1, m)
        subaxis(size(section_sort, 1), 1, section_id, 'Spacing', 0.03, 'Padding', 0, 'Margin', 0);
        
        section_now = [section_sort(section_id, 1), section_sort(section_id, 2)];
        section_flag = spike_time >= section_now(1) & spike_time <= section_now(2);

        spike_time_section = spike_time(section_flag);
        spike_time_section = spike_time_section - section_now(1);

        rep_len = 3;
        rep_max = round(section_now(2) - section_now(1)) / rep_len;

        for rep = 1 : rep_max
            rep_flag = spike_time_section >= (rep-1)*rep_len & spike_time_section <= rep*rep_len;
            spike_time_rep = spike_time_section(rep_flag);
            spike_time_rep = spike_time_rep - (rep-1)*rep_len;

            rep_mark = rep * ones(length(spike_time_rep),1);
            scatter(spike_time_rep, rep_mark, 5, 'filled')
            axis([-0.05 (rep_len + 0.05) 0 (rep_max + 1)])
            hold on
        end
        hold on
        axis off
    end
        
    saveas(gcf, [num2str(ds_slave_id), '-unsorted-pca-stitch.jpg'])
    savefig([num2str(ds_slave_id), '-unsorted-pca-stitch.fig'])
    disp(['saved fig for ', num2str(ds_slave_id)])
    close
end
    
toc % takes 7-9 min to generate a single cell sensitivity plot. needs optim
