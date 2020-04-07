%% load data 

clear
clc

dataset_num = '/data000-map'; 
date_num = '2016-03-04-0';

% prefix_now = '/Volumes/dusom_fieldlab';
prefix_now = '/Volumes/All_Staff/';

datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, dataset_num, dataset_num);
datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);
% datarun = load_ei(datarun, 'all', 'array_type', 519);

ds_slave_id_seq = importdata('slave_2016.txt'); % 20160304 data

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
nflash = [150; sections(2:end,3)./flash_period];
sections = [sections, ndf, flash_config, nflash];

sections(:,7) = (sections(:,4)*(-10) + sections(:,5));
section_sort = sortrows(sections, 7);

binsize = 0.020; 
trial_len = 3 / binsize; % all trial length = 3s
binnum = datarun.duration / binsize;
edges = linspace(0, datarun.duration, binnum);
section_idx = [round(section_sort(:,1)/binsize), round(section_sort(:,2)/binsize), section_sort];
section_idx(:,end) = (section_idx(:,6) * 10 + section_idx(:,7));
ntrial = round(section_idx(:,8));
marker = unique(section_idx(:,end), 'stable');

%% convert x axis to Rh*/rod

x0 = 854; xunit = -251; r = 13;
right_edge = [40,120,196,301, 344,376,420,452, 544,620,694,774,927];
xpos = right_edge - r;
x = (x0 - xpos)' / xunit;
intensity = 10.^x;

% scatter(x, ones(length(x),1))
% xline(-3); xline(-2); xline(-1); xline(0);
% conjecture: used marker(2:end-2), aka from 53.2 to 23.4

x_n_marker = [x, marker(2:14)];
% save xm x_n_marker

%% iterate across cells

cell_excluded = 0;
for c = 1 : 2 %length(ds_slave_id_seq)
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
        for flash_intensity = 2 : length(marker) % exclude dark==990
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
%         x = 1 : length(marker)-1;
        errorbar(x, Pc_avg(1:end-2), Pc_var(1:end-2))
        hold on
        yline(1,'-.g'); yline(0.84,'-.g');
%         xticks(x)
%         xticklabels({'5.002','5.004','5.008','4.002','4.003','4.004','4.006','4.008',...
%             '3.002','3.004','3.008','2.002','2.004','1.002','0.002'})
%         xtickangle(45)
        xlabel('log(intensity)')
        ylabel('probability correct')
        ylim([0.45, 1.03])
        
%         saveas(gcf, ['log_intensity-', num2str(ds_slave_id_seq(c)), '.png'])
%         print(['log_intensity-', num2str(ds_slave_id_seq(c))], '-dpdf', '-fillpage')
%         disp(['saved fig for ', num2str(ds_slave_id_seq(c))])
%         close
    else
        cell_excluded = cell_excluded + 1;
        disp(['cell excluded ', num2str(ds_slave_id_seq(c))])
        close
    end
end
