% sensitivity plot for single cell
% later will import ds cell id list to achieve sequential plotting for
% every single ds cell

%% load data w changeable dataset_num & ds_now
clear
clc

dataset_num = '00'; % dim flashes stimulus to test absolute sensitivity
prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/2019-11-21-0/data0', dataset_num, '/data0', dataset_num);
% datapath = '/Volumes/???/lab/Experiments/Array/Analysis/2019-11-21-0/data002/data002';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
datarun = load_ei(datarun, 'all', 'array_type', 519);

%% chop data000 to sections
gaps = round(diff(datarun.triggers));
switch_flag = gaps ~= 2 & gaps ~= 4;
switch_duration = [round(datarun.triggers(1)) - 300; gaps(switch_flag)];

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

%% scaffold: single section of darkness. cell id = 469, index = 21
cell_index = 35; % slave index! find by datarun.cell_ids. match master id to slave id to slave index
spike_time = datarun.spikes{cell_index, 1};

section_id = 10; % range 1-17
section_now = [sections(section_id, 1), sections(section_id, 2)];

section_flag = spike_time >= section_now(1) & spike_time <= section_now(2);
spike_time_section = spike_time(section_flag);
spike_time_section = spike_time_section - section_now(1);

rep_max = round(section_now(2) - section_now(1)) / 2;
for rep = 1 : rep_max
    rep_flag = spike_time_section >= (rep-1)*2 & spike_time_section <= rep*2;
    spike_time_rep = spike_time_section(rep_flag);
    spike_time_rep = spike_time_rep - (rep-1)*2;

    rep_mark = rep * ones(length(spike_time_rep),1);
    scatter(spike_time_rep, rep_mark, 10, 'filled')
    hold on
    axis([-0.05 2.05 0 (rep_max + 1)])
end

%% plot all sections separately
for section_i = 1 : size(sections, 1)
    subplot( size(sections, 1), 1, section_i )
    title(['NDF = ', num2str(sections(section_i, 4))])

    section_id = section_i; 
    section_now = [sections(section_id, 1), sections(section_id, 2)];

    section_flag = spike_time >= section_now(1) & spike_time <= section_now(2);
    spike_time_section = spike_time(section_flag);
    spike_time_section = spike_time_section - section_now(1);

    rep_max = round(section_now(2) - section_now(1)) / 2;
    for rep = 1 : rep_max
        rep_flag = spike_time_section >= (rep-1)*2 & spike_time_section <= rep*2;
        spike_time_rep = spike_time_section(rep_flag);
        spike_time_rep = spike_time_rep - (rep-1)*2;

        rep_mark = rep * ones(length(spike_time_rep),1);
        scatter(spike_time_rep, rep_mark, 10, 'filled')
        hold on
        axis([-0.05 2.05 0 (rep_max + 1)])
    end
end

%% sort sections by NDF and flash_config. fixed x_axis=2
sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);

for section_i = 1 : size(section_sort, 1)
    subplot( size(section_sort, 1), 1, section_i )

    section_id = section_i; 
    section_now = [section_sort(section_id, 1), section_sort(section_id, 2)];

    section_flag = spike_time >= section_now(1) & spike_time <= section_now(2);
    spike_time_section = spike_time(section_flag);
    spike_time_section = spike_time_section - section_now(1);
    
    if section_i >=1 && section_i <= 15
        rep_len = 2;
    else
        rep_len = 4;
    end

    rep_max = round(section_now(2) - section_now(1)) / rep_len;
    for rep = 1 : rep_max
        rep_flag = spike_time_section >= (rep-1)*rep_len & spike_time_section <= rep*rep_len;
        spike_time_rep = spike_time_section(rep_flag);
        spike_time_rep = spike_time_rep - (rep-1)*rep_len;

        rep_mark = rep * ones(length(spike_time_rep),1);
        scatter(spike_time_rep, rep_mark, 10, 'filled')
        hold on
        axis([-0.05 2.05 0 (rep_max + 1)])
    end
%     title(['NDF = ', num2str(section_sort(section_i, 4))])
end

%% merge sections w same NDF and flash_config. flexible x_axis=2 or 4
marker = unique(section_sort(:,7));
marker_seq = section_sort(:,7);

for m = 1:length(marker)
    subplot(length(marker), 1, m)  
    marker_now = marker(m);
    section_id_seq = find(marker_seq == marker_now, length(marker_seq));
    
    for s = 1:length(section_id_seq)
        section_id = section_id_seq(s); 
        section_now = [section_sort(section_id, 1), section_sort(section_id, 2)];
        section_flag = spike_time >= section_now(1) & spike_time <= section_now(2);
        
        spike_time_section = spike_time(section_flag);
        spike_time_section = spike_time_section - section_now(1);

        if section_i >=1 && section_i <= 15
            rep_len = 2;
        else
            rep_len = 2;
%             rep_len = 4;
        end
        rep_max = round(section_now(2) - section_now(1)) / rep_len;
        
        for rep = 1 : rep_max
            rep_flag = spike_time_section >= (rep-1)*rep_len & spike_time_section <= rep*rep_len;
            spike_time_rep = spike_time_section(rep_flag);
            spike_time_rep = spike_time_rep - (rep-1)*rep_len;

            rep_mark = rep * ones(length(spike_time_rep),1);
            scatter(spike_time_rep, rep_mark, 10, 'filled')
            hold on
            axis([-0.05 (rep_len + 0.05) 0 (rep_max + 1)])
        end
        hold on
    end    
end
