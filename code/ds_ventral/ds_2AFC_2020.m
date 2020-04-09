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

ds_master_slave_id = importdata('20200229_ds_map_final.txt'); 
ds_slave_id_seq = ds_master_slave_id(:,2); % 20200229 data

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

binsize = 0.020; 
trial_len = 2 / binsize; % pretend all trial length = 2s
binnum = datarun.duration / binsize;
edges = linspace(0, datarun.duration, binnum);
section_idx = [round(section_sort(:,1)/binsize), round(section_sort(:,2)/binsize), section_sort];
section_idx(:,end) = (section_idx(:,6) * 10 + section_idx(:,7));
marker = unique(section_idx(:,end), 'stable');
ntrial = round(section_idx(:,8));

%%
load('xm.mat')
marker_match = marker(2:end) - mod(floor(marker(2:end)),10) + 3;
flash_1ms = abs(marker_match - floor(marker_match))-0.1 < 1e-4;
marker_match = marker_match(~flash_1ms);

temp_id = ismembertol(x_n_marker(:,2), marker_match, 1e-4);
x = x_n_marker(temp_id,1);

%% iterate across cells

cell_excluded = 0;
for c = 1 : 2 %1 : length(ds_slave_id_seq)
    figure
    ds_slave_index = find(datarun.cell_ids == ds_slave_id_seq(c)); 
    spike_time = datarun.spikes{ds_slave_index, 1};
    [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin
    
    sum_null = zeros(ntrial(1), trial_len); 
    for t = 1 : ntrial(1)
        trial_null = binned(trial_len*(t-1)+section_idx(1,1)+1 : trial_len*t+section_idx(1,1));
        sum_null(t,:) = trial_null;
    end
    sum_null = sum(sum_null,1);

    ntest = 1000;
    Pc = zeros(length(marker)-1 ,ntest);
    for test = 1 : ntest
        for flash_intensity = 2 : length(marker) % exclude dark==990. should improve by excluding 1ms here
            nid = 1; 
            fid = section_idx(:,end)==marker(flash_intensity);
            fid_seq = find(fid==1);
            
            if floor(section_idx(fid_seq(1),7)) == 4
                scale = 2; % account for 4s trials
            else
                scale = 1;
            end

            sum_flash_seq = zeros(1, trial_len*scale);
            for i = 1 : length(fid_seq)
                sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len*scale);
                for t = 1 : ntrial(fid_seq(i))
                    trial_flash = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
                    sum_flash_section{i}(t,:) = trial_flash;
                end
                sum_flash_section{i} = sum(sum_flash_section{i},1);
                sum_flash_seq = sum_flash_seq + sum_flash_section{i};
            end
            sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
            
            sum_all = sum_flash_seq + sum_null;
            trial_num_null = 1 : ntrial(nid);
            trial_num_flash = 1 : scale : scale*sum(ntrial(fid));
            mean_all = sum_all ./ (length(trial_num_null) + length(trial_num_null));
            
            sample_size = min(length(trial_num_null), length(trial_num_flash));
            order_null = datasample(trial_num_null, sample_size, 'Replace', false);
            order_flash = datasample(trial_num_flash, sample_size, 'Replace', false);

            corrpos = zeros(sample_size, 1);
            for t = 1 : sample_size
                trial_null = binned(trial_len*(order_null(t)-1)+section_idx(nid,1)+1 : trial_len*order_null(t)+section_idx(nid,1));
                other_null = sum_null - trial_null;
                mean_null = other_null ./ (length(trial_num_null) - 1) - mean_all .* length(trial_num_null) ./ (length(trial_num_null) - 1); % zero-mean

                if scale == 1
                    if order_flash(t) <= ntrial(fid_seq(1))
                        trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_flash(t)+section_idx(fid_seq(1),1));
                    else
                        trial_flash = binned(trial_len*(order_flash(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                            trial_len*(order_flash(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
                    end
                elseif scale == 2
                    trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_flash(t)+section_idx(fid_seq(1),1));
                end
                other_flash = sum_flash_seq - trial_flash;
                mean_flash = other_flash ./ (length(trial_num_flash) - 1) - mean_all .* length(trial_num_flash) ./ (length(trial_num_flash) - 1);

                discriminant = (mean_flash - mean_null)';
                corrpos(t) = (trial_flash - trial_null) * discriminant; % no need to zero mean trial_flash & _null because they cancel out
            end
            corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
            Pc(flash_intensity - 1, test) = corr / length(corrpos);
        end
    end

    Pc_avg = mean(Pc,2);
    Pc_var = std(Pc,1,2);
    
    if max(Pc_avg) >= 0.84
%         x = 1 : length(marker)-3; % exclude dark & 1ms flash (990, 52.1, 42.1)
%         flash_1ms = abs(marker(2:end) - floor(marker(2:end)))-0.1 < 1e-4;
        Pc_avg = Pc_avg(~flash_1ms);
        Pc_var = Pc_var(~flash_1ms);

%         x = 1 : length(marker)-1;
        errorbar(x, Pc_avg, Pc_var)
        hold on
        yline(1,'-.g'); yline(0.84,'-.g');
%         xticks(x)
%         xticklabels({'52.1','52.2','52.4','52.8','42.1','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})
%         xticklabels({'52.2','52.4','52.8','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})        
%         xtickangle(45)
        xlabel('log(intensity)')
        ylabel('probability correct')
        ylim([0.4, 1.05])
        
%         saveas(gcf, ['log_intensity-', num2str(ds_slave_id_seq(c)), '.png'])
%         print(['log_intensity-', num2str(ds_slave_id_seq(c))], '-dpdf', '-fillpage')
%         disp(['saved fig for ', num2str(ds_slave_id_seq(c))])
%         close
    else
        cell_excluded = cell_excluded + 1;
        disp([num2str(ds_slave_id_seq(c)),' excluded due to low Pc'])
%         close
    end
end

%% test Discriminant shape

for c = 1 
%     figure
    ds_slave_index = find(datarun.cell_ids == ds_slave_id_seq(c)); 
    spike_time = datarun.spikes{ds_slave_index, 1};
    [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin

    sum_null = zeros(ntrial(1), trial_len); 
    for t = 1 : ntrial(1)
        trial_null = binned(trial_len*(t-1)+section_idx(1,1)+1 : trial_len*t+section_idx(1,1));
        sum_null(t,:) = trial_null;
    end
    sum_null = sum(sum_null,1);
    D_all = cell(length(marker),1);

    ntest = 1000;
    Pc = zeros(length(marker)-1 ,ntest);
    for test = 1 : ntest
        for flash_intensity = 2 : length(marker) % exclude dark==990. should improve by excluding 1ms here
            nid = 1; 
            fid = section_idx(:,end)==marker(flash_intensity);
            fid_seq = find(fid==1);
            D_fid = zeros(ntest*120, 100);
            
            if floor(section_idx(fid_seq(1),7)) == 4
                scale = 2; % account for 4s trials
            else
                scale = 1;
            end

            sum_flash_seq = zeros(1, trial_len*scale);
            for i = 1 : length(fid_seq)
                sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len*scale);
                for t = 1 : ntrial(fid_seq(i))
                    trial_flash = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
                    sum_flash_section{i}(t,:) = trial_flash;
                end
                sum_flash_section{i} = sum(sum_flash_section{i},1);
                sum_flash_seq = sum_flash_seq + sum_flash_section{i};
            end
            sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
            
            sum_all = sum_flash_seq + sum_null;
            trial_num_null = 1 : ntrial(nid);
            trial_num_flash = 1 : scale : scale*sum(ntrial(fid));
            mean_all = sum_all ./ (length(trial_num_null) + length(trial_num_null));
            
            sample_size = min(length(trial_num_null), length(trial_num_flash));
            order_null = datasample(trial_num_null, sample_size, 'Replace', false);
            order_flash = datasample(trial_num_flash, sample_size, 'Replace', false);

            corrpos = zeros(sample_size, 1);
            for t = 1 : sample_size
                trial_null = binned(trial_len*(order_null(t)-1)+section_idx(nid,1)+1 : trial_len*order_null(t)+section_idx(nid,1));
                other_null = sum_null - trial_null;
                mean_null = other_null ./ (length(trial_num_null) - 1) - mean_all .* length(trial_num_null) ./ (length(trial_num_null) - 1); % zero-mean

                if scale == 1
                    if order_flash(t) <= ntrial(fid_seq(1))
                        trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_flash(t)+section_idx(fid_seq(1),1));
                    else
                        trial_flash = binned(trial_len*(order_flash(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                            trial_len*(order_flash(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
                    end
                elseif scale == 2
                    trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_flash(t)+section_idx(fid_seq(1),1));
                end
                other_flash = sum_flash_seq - trial_flash;
                mean_flash = other_flash ./ (length(trial_num_flash) - 1) - mean_all .* length(trial_num_flash) ./ (length(trial_num_flash) - 1);

                discriminant = (mean_flash - mean_null)';
                D_fid((test-1)*60 + t, :) = discriminant;
                corrpos(t) = (trial_flash - trial_null) * discriminant; % no need to zero mean trial_flash & _null because they cancel out
            end
            corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
            Pc(flash_intensity - 1, test) = corr / length(corrpos);
            D_all{flash_intensity,1} = mean(D_fid,1);

        end
    end


    Pc_avg = mean(Pc,2);
    Pc_var = std(Pc,1,2);
    
    if max(Pc_avg) >= 0.84

        Pc_avg = Pc_avg(~flash_1ms);
        Pc_var = Pc_var(~flash_1ms);

%         errorbar(x, Pc_avg, Pc_var)
%         hold on
%         yline(1,'-.g'); yline(0.84,'-.g');
%         xticks(x)
%         xticklabels({'52.2','52.4','52.8','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})        
%         xtickangle(45)
% %         xlabel('log(intensity)')
% %         ylabel('probability correct')
%         ylim([0.4, 1.05])
    
    end
end

%%
figure('units','normalized','outerposition',[0 0 1 1]) 
for i = 2:13
%     subplot(12,1,i-1)
    subaxis(12,1,i-1, 'Spacing', 0.001, 'Padding', 0, 'Margin', 0);
    D_intensity = D_all{i,1};
    plot(D_intensity)
end