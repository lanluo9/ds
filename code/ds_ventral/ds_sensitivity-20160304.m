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

% %% load ds cell identified in master 

% % load('ds_cell_map_20200306.mat', 'ds_map_all');
% flag = find(ds_map_all(:,1)==0);
% ds_slave_id_mapPCA = ds_map_all(1:(flag(1)-1), 2); ismember(ds_slave_id_mapPCA, datarun.cell_ids)
% ds_slave_id_mapEI = ds_map_all((flag(1)+1):(flag(2)-1), 2);
% ds_slave_id_map2 = ds_map_all((flag(2)+1):end, 2); ds_slave_id_map2(ds_slave_id_map2 == 0) = [];
% % ds_map_all

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
sections(end+1, :) = [datarun.triggers(end)+flash_period, datarun.duration, datarun.duration - (datarun.triggers(end)+flash_period)];

ndf = [99; 5;5;5; 4;4;4;4;4; 3;3;3; 2;2; 1;0; 99];
flash_config = [0; 3.2;3.4;3.8; 3.2;3.3;3.4;3.6;3.8; 3.2;3.4;3.8; 3.2;3.4; 3.2;3.2; 0];
nflash = [0; sections(2:end-1,3)./flash_period; 0];
sections = [sections, ndf, flash_config, nflash];

sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);
marker = unique(section_sort(:,7));
marker_seq = section_sort(:,7);

%% merge sections w same NDF and flash_config. x_axis=2 after cutting off 2-4s

slave_ds_id_all = unique(ds_map_all(:,2)); slave_ds_id_all(slave_ds_id_all == 0) = [];

% slave_ds_id_all = [3318 3664 3707 6841 3061 5629 5495]; % lost bc of slave cleaning
% 
% for i = 1 : length(slave_ds_id_all)
%     ds_slave_id = slave_ds_id_all(i); 
%     ds_slave_index = find(datarun.cell_ids == ds_slave_id); 
%     if ~isempty(ds_slave_index)
%         disp([num2str(ds_slave_id), '  found in slave datarun.cell_id'])
%         continue
%     end
% end

tic

for i = 1 : length(slave_ds_id_all)
    figure 

    ds_slave_id = slave_ds_id_all(i); 
    ds_slave_index = find(datarun.cell_ids == ds_slave_id); 
    if isempty(ds_slave_index)
        disp([num2str(ds_slave_id), ' not found in slave datarun.cell_id'])
        continue
    end
    
    spike_time = datarun.spikes{ds_slave_index, 1};

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
    
%     title(['data0', num2str(dataset_num), '. dsRGC index = ', num2str(ds_slave_index), '. id = ', num2str(ds_slave_id)])
    saveas(gcf, ['map_lost-', num2str(ds_slave_index),'-', num2str(ds_slave_id), '.png'])
    close
end
    
toc % takes 7-9 min to generate a single cell sensitivity plot. needs optim
