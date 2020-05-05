%% load data 

clear
clc
close

% dataset_num = '/data000-map-sorted'; 
% date_num = '2020-02-29-0';
% % prefix_now = '/Volumes/dusom_fieldlab';
% prefix_now = '/Volumes/All_Staff/';
% datapath = append(prefix_now, '/lab/Experiments/Array/Analysis/', date_num, dataset_num, dataset_num);

datapath = 'D:/RRR/Grad/Rotation/GF_lab/lab_Mac/2020-02-29-0/data000-map-sorted/data000-map-sorted';

datarun = load_data(datapath);
datarun = load_neurons(datarun);
datarun = load_params(datarun);

% ds_master_slave_id = importdata('20200229_ds_map_final.txt'); 
% ds_slave_id_seq = ds_master_slave_id(:,2); % 20200229 data

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

% %%
% load('xm.mat')
% marker_match = marker(2:end) - mod(floor(marker(2:end)),10) + 3;
% flash_1ms = abs(marker_match - floor(marker_match))-0.1 < 1e-4;
% marker_match = marker_match(~flash_1ms);
% 
% temp_id = ismembertol(x_n_marker(:,2), marker_match, 1e-4);
% x = x_n_marker(temp_id,1);

%% load latency txt & x coor
normal_latency = importdata('latency.txt')
divide = find(normal_latency(:,1)==0);
ds_slave_normal = normal_latency(1:(divide-1), 2);
ds_slave_latency = normal_latency((divide+1):end, 2);
ds_slave_all = [ds_slave_normal; ds_slave_latency];

load('xm.mat')
marker_match = marker(2:end) - mod(floor(marker(2:end)),10) + 3;
flash_1ms = abs(marker_match - floor(marker_match))-0.1 < 1e-4;
marker_match = marker_match(~flash_1ms);
temp_id = ismembertol(x_n_marker(:,2), marker_match, 1e-4);
log_intensity = x_n_marker(temp_id,1);

%% test simple spike counter (0-300 ms after flash)

binsize = 0.020;
trial_len = 2 / binsize; % pretend all trial length = 2s
binnum = datarun.duration / binsize;
edges = linspace(0, datarun.duration, binnum);
section_idx = [round(section_sort(:,1)/binsize), round(section_sort(:,2)/binsize), section_sort];
section_idx(:,end) = (section_idx(:,6) * 10 + section_idx(:,7));
marker = unique(section_idx(:,end), 'stable');
ntrial = round(section_idx(:,8));

counter_range = [0.200,0.300,0.400]; % expect ~200 to work
pc = struct;

ds_slave_now = sort(ds_slave_normal);

for r = 1 : length(counter_range)    
    range = counter_range(r);
    
    for cellnum = 1 : length(ds_slave_now)
        
        n = (r-1)*length(ds_slave_now) + cellnum;
        pc(n).range = range;
        pc(n).cellid = ds_slave_now(cellnum);
        disp(['calculating for cell ' num2str(ds_slave_now(cellnum)) ' counter_range ' num2str(range)])
        
        ds_slave_index = find(datarun.cell_ids == ds_slave_now(cellnum)); 
        spike_time = datarun.spikes{ds_slave_index, 1};
        [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin

        ntest = 10;
        Pc = zeros(length(marker)-1 ,ntest);
        for test = 1 : ntest
            for flash_intensity = length(marker) :-1: 2 % exclude dark==990. should improve by excluding 1ms here
                fid = section_idx(:,end)==marker(flash_intensity);
                fid_seq = find(fid==1);

                if floor(section_idx(fid_seq(1),7)) == 4, scale = 2; % account for 4s trials
                else, scale = 1; end

%                 sum_flash_seq = zeros(1, trial_len*scale);
%                 for i = 1 : length(fid_seq)
%                     sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len*scale);
%                     for t = 1 : ntrial(fid_seq(i))
%                         trial_flash_full = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
%                         sum_flash_section{i}(t,:) = trial_flash_full;
%                     end
%                     sum_flash_section{i} = sum(sum_flash_section{i},1);
%                     sum_flash_seq = sum_flash_seq + sum_flash_section{i};
%                 end
%                 sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
%                 sum_flash_pre = sum_flash_seq(1 : length(sum_flash_seq)/2); % pre (0-1s) as flash trial
%                 sum_flash_post = sum_flash_seq(length(sum_flash_seq)/2 + 1 : end); % post (1-2s) as null trial
% 
%                 sum_all = sum_flash_pre + sum_flash_post;
                trial_num_post = 1 : scale : scale*sum(ntrial(fid));
                trial_num_pre = 1 : scale : scale*sum(ntrial(fid));
%                 mean_all = sum_all ./ (length(trial_num_post) + length(trial_num_pre)); 

                sample_size = min(length(trial_num_post), length(trial_num_pre));
                order_post = datasample(trial_num_post, sample_size, 'Replace', false);
                order_pre = datasample(trial_num_pre, sample_size, 'Replace', false);

                corrpos = zeros(sample_size, 1);
                for t = 1 : sample_size
                    if scale == 1
                        if order_pre(t) <= ntrial(fid_seq(1))
                            trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                                trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1) + range/binsize);
                        else
                            trial_flash_pre = binned(trial_len*(order_pre(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                                trial_len*(order_pre(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1) + range/binsize);
                        end
                    elseif scale == 2
                        trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                                trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1) + range/binsize);
                    end
                    counter_pre = sum(trial_flash_pre(:));
%                     other_flash_pre = sum_flash_pre - trial_flash_pre;
%                     mean_flash_pre = other_flash_pre ./ (length(trial_num_pre) - 1) - mean_all .* length(trial_num_pre) ./ (length(trial_num_pre) - 1);

                    if scale == 1
                        if order_post(t) <= ntrial(fid_seq(1))
                            trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 :...
                                trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + range/binsize);
                        else
                            trial_flash_post = binned(trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2 + 1 : ...
                                trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2 + range/binsize);
                        end
                    elseif scale == 2 
                        trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 : ...
                                trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + range/binsize);
                    end
                    counter_post = sum(trial_flash_post(:));
%                     other_flash_post = sum_flash_post - trial_flash_post;
%                     mean_flash_post = other_flash_post ./ (length(trial_num_post) - 1) - mean_all .* length(trial_num_post) ./ (length(trial_num_post) - 1);

%                     discriminant = (mean_flash_pre - mean_flash_post)';
                    corrpos(t) = counter_pre - counter_post;
                end
                corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
                Pc(flash_intensity - 1, test) = corr / length(corrpos);
            end
        end

        Pc_avg = mean(Pc,2);
        Pc_var = std(Pc,1,2);
        Pc_avg = Pc_avg(~flash_1ms)
        Pc_var = Pc_var(~flash_1ms);
        pc(n).Pc_avg = Pc_avg;
        pc(n).Pc_var = Pc_var;
    end
    disp(['done calculating for counter_range ' num2str(range)])
end
disp('calc done')

%% counter range impact on Pc

T = struct2table(pc);
sortedT = sortrows(T, 'cellid');
pc_sorted = table2struct(sortedT); % dependent on: ds_slave_now = sort(ds_slave_now);
% save pc_sorted

color = prism(length(counter_range));
disp('start plotting')

for cellnum = 1 : length(ds_slave_now)
    for r = 1 : length(counter_range)
        n = (cellnum - 1) * length(counter_range) + r;
        Pc_avg = pc_sorted(n).Pc_avg;
%         Pc_var = pc_sorted(n).Pc_var;
        plot(log_intensity, Pc_avg, 'Color', color(r,:), 'LineWidth', 2)
%         errorbar(log_intensity, Pc_avg, Pc_var, 'Color', color(b,:), 'LineWidth', 2)
        hold on
    end
    line([min(log_intensity), max(log_intensity)], [0.84, 0.84],'Color', [0 1 0])
    
    ylim([0.4, 1.05])
    xlabel('log(intensity)')
    ylabel('probability correct')
    legend({'200','300','400'},'Location','southeast')
    legend('boxoff')
    
    saveas(gcf, [num2str(ds_slave_now(cellnum)) '-2AFC-fit-counter-range' '.png'])
    disp(['saved fig for ', num2str(ds_slave_now(cellnum))])
    close
end    

%% test binsize

binsize_seq = [0.020, 0.050, 0.100, 0.200, 0.250, 0.500]; % 0.300 yields non-int trial_len
pc = struct;

ds_slave_now = ds_slave_normal;

for b = 1 : length(binsize_seq)
    
    binsize = binsize_seq(b); % bin window = 20 ms - 500 ms
    
    trial_len = 2 / binsize; % pretend all trial length = 2s
    binnum = datarun.duration / binsize;
    edges = linspace(0, datarun.duration, binnum);
    section_idx = [round(section_sort(:,1)/binsize), round(section_sort(:,2)/binsize), section_sort];
    section_idx(:,end) = (section_idx(:,6) * 10 + section_idx(:,7));
    marker = unique(section_idx(:,end), 'stable');
    ntrial = round(section_idx(:,8));

    for cellnum = 1 : length(ds_slave_now)
        
        n = (b-1)*length(ds_slave_now) + cellnum;
        pc(n).binsize = binsize;
        pc(n).cellid = ds_slave_now(cellnum);
        disp(['calculating for cell ' num2str(ds_slave_now(cellnum)) ' binsize ' num2str(binsize)])
        
    %     figure('units','normalized','outerposition',[0 0 1 1]) 
        ds_slave_index = find(datarun.cell_ids == ds_slave_now(cellnum)); 
        spike_time = datarun.spikes{ds_slave_index, 1};
        [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin

        ntest = 1000;
        Pc = zeros(length(marker)-1 ,ntest);
        for test = 1 : ntest
            for flash_intensity = 2 : length(marker) % exclude dark==990. should improve by excluding 1ms here
                fid = section_idx(:,end)==marker(flash_intensity);
                fid_seq = find(fid==1);

                if floor(section_idx(fid_seq(1),7)) == 4, scale = 2; % account for 4s trials
                else, scale = 1; end

                sum_flash_seq = zeros(1, trial_len*scale);
                for i = 1 : length(fid_seq)
                    sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len*scale);
                    for t = 1 : ntrial(fid_seq(i))
                        trial_flash_full = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
                        sum_flash_section{i}(t,:) = trial_flash_full;
                    end
                    sum_flash_section{i} = sum(sum_flash_section{i},1);
                    sum_flash_seq = sum_flash_seq + sum_flash_section{i};
                end
                sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
                sum_flash_pre = sum_flash_seq(1 : length(sum_flash_seq)/2); % pre (0-1s) as flash trial
                sum_flash_post = sum_flash_seq(length(sum_flash_seq)/2 + 1 : end); % post (1-2s) as null trial

                sum_all = sum_flash_pre + sum_flash_post;
                trial_num_post = 1 : scale : scale*sum(ntrial(fid));
                trial_num_pre = 1 : scale : scale*sum(ntrial(fid));
                mean_all = sum_all ./ (length(trial_num_post) + length(trial_num_pre)); 

                sample_size = min(length(trial_num_post), length(trial_num_pre));
                order_post = datasample(trial_num_post, sample_size, 'Replace', false);
                order_pre = datasample(trial_num_pre, sample_size, 'Replace', false);

                corrpos = zeros(sample_size, 1);
                for t = 1 : sample_size
                    if scale == 1
                        if order_pre(t) <= ntrial(fid_seq(1))
                            trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                                trial_len*order_pre(t)+section_idx(fid_seq(1),1) - trial_len/2);
                        else
                            trial_flash_pre = binned(trial_len*(order_pre(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                                trial_len*(order_pre(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2);
                        end
                    elseif scale == 2
                        trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                                trial_len*order_pre(t)+section_idx(fid_seq(1),1) - trial_len/2);
                    end
                    other_flash_pre = sum_flash_pre - trial_flash_pre;
                    mean_flash_pre = other_flash_pre ./ (length(trial_num_pre) - 1) - mean_all .* length(trial_num_pre) ./ (length(trial_num_pre) - 1);

                    if scale == 1
                        if order_post(t) <= ntrial(fid_seq(1))
                            trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 :...
                                trial_len*order_post(t)+section_idx(fid_seq(1),1));
                        else
                            trial_flash_post = binned(trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2 + 1 : ...
                                trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
                        end
                    elseif scale == 2 
                        trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 : ...
                                trial_len*order_post(t)+section_idx(fid_seq(1),1));
                    end
                    other_flash_post = sum_flash_post - trial_flash_post;
                    mean_flash_post = other_flash_post ./ (length(trial_num_post) - 1) - mean_all .* length(trial_num_post) ./ (length(trial_num_post) - 1);

                    discriminant = (mean_flash_pre - mean_flash_post)';
                    corrpos(t) = (trial_flash_pre - trial_flash_post) * discriminant; % no need to zero mean trial_flash & _null because they cancel out
                end
                corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
                Pc(flash_intensity - 1, test) = corr / length(corrpos);
            end
        end

        Pc_avg = mean(Pc,2);
        Pc_var = std(Pc,1,2);
        Pc_avg = Pc_avg(~flash_1ms);
        Pc_var = Pc_var(~flash_1ms);
        pc(n).Pc_avg = Pc_avg;
        pc(n).Pc_var = Pc_var;
%         c = x;
%         r = Pc_avg;
%         fit = fitNakaRushton(c,r)

%         hold on
%         c = linspace(min(x), max(x), 10^3);
%         Pc_fit = fit.Rmax * ((c.^fit.n) ./ ((c.^fit.n) + fit.c50.^fit.n)) + fit.offset;
%         line([min(x), max(x)], [1, 1],'Color', [0 1 0])
%         line([min(x), max(x)], [0.84, 0.84],'Color', [0 1 0])
%         plot(c, Pc_fit, 'b--', 'LineWidth', 1)
%         errorbar(x, Pc_avg, Pc_var, 'Color', [1 0 0])
% 
%         ylim([0.4, 1.05])
%         titleStr = sprintf('c50: %0.2f n: %0.2f\n  Rmax: %0.2f offset: %0.2f', fit.c50,fit.n, fit.Rmax,fit.offset);
%         title(sprintf(titleStr));
%         xlabel('log(intensity)')
%         ylabel('probability correct')
%         ylim([0.4, 1.05])
%         saveas(gcf, [num2str(ds_slave_now(cellnum)) '-2AFC-fit-' num2str(1000*binsize_seq(b)) '.png'])
%         disp(['saved fig for ', num2str(ds_slave_now(cellnum))])
%         close
    end
    disp(['done calculating for binsize ' num2str(binsize)])
end
disp('calc done')

%% binsize impact on Pc

T = struct2table(pc);
sortedT = sortrows(T, 'cellid');
pc_sorted = table2struct(sortedT);
save pc_sorted

color = prism(6);
disp('start plotting')

ds_slave_now = sort(ds_slave_now);
for cellnum = 1 : length(ds_slave_now)
    for b = 1 : length(binsize_seq)
        n = (cellnum - 1) * length(binsize_seq) + b;
        Pc_avg = pc_sorted(n).Pc_avg;
        Pc_var = pc_sorted(n).Pc_var;
        plot(log_intensity, Pc_avg, 'Color', color(b,:), 'LineWidth', 2)
%         errorbar(log_intensity, Pc_avg, Pc_var, 'Color', color(b,:), 'LineWidth', 2)
        hold on
    end
    line([min(log_intensity), max(log_intensity)], [0.84, 0.84],'Color', [0 1 0])
    
    ylim([0.4, 1.05])
    xlabel('log(intensity)')
    ylabel('probability correct')
    legend({'20','50','100','200','250','500'},'Location','southeast')
    legend('boxoff')
    
    saveas(gcf, [num2str(ds_slave_now(cellnum)) '-2AFC-fit-binsize' '.png'])
    disp(['saved fig for ', num2str(ds_slave_now(cellnum))])
    close
end    


%% use 1-2s as null trial for normal latency cells

ds_slave_all = sort(ds_slave_all);
ds_slave_now = ds_slave_all;
% cell_excluded = 0;
for cellnum = 1 : length(ds_slave_now)
%     figure('units','normalized','outerposition',[0 0 1 1]) 
    
    if ismember(ds_slave_now(cellnum), ds_slave_normal), flag = 'normal';
    elseif ismember(ds_slave_now(cellnum), ds_slave_latency), flag = 'long latency'; end

    ds_slave_index = find(datarun.cell_ids == ds_slave_now(cellnum)); 
    spike_time = datarun.spikes{ds_slave_index, 1};
    [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin
    
    ntest = 1000;
    Pc = zeros(length(marker)-1 ,ntest);
    for test = 1 : ntest
        for flash_intensity = 2 : length(marker) % exclude dark==990. should improve by excluding 1ms here
%             nid = 1; 
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
                    trial_flash_full = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
                    sum_flash_section{i}(t,:) = trial_flash_full;
                end
                sum_flash_section{i} = sum(sum_flash_section{i},1);
                sum_flash_seq = sum_flash_seq + sum_flash_section{i};
            end
            sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
            sum_flash_pre = sum_flash_seq(1 : length(sum_flash_seq)/2); % pre (0-1s) as flash trial
            sum_flash_post = sum_flash_seq(length(sum_flash_seq)/2 + 1 : end); % post (1-2s) as null trial
            
            sum_all = sum_flash_pre + sum_flash_post;
            trial_num_post = 1 : scale : scale*sum(ntrial(fid));
            trial_num_pre = 1 : scale : scale*sum(ntrial(fid));
            mean_all = sum_all ./ (length(trial_num_post) + length(trial_num_pre)); 
            
            sample_size = min(length(trial_num_post), length(trial_num_pre));
            order_post = datasample(trial_num_post, sample_size, 'Replace', false);
            order_pre = datasample(trial_num_pre, sample_size, 'Replace', false);

            corrpos = zeros(sample_size, 1);
            for t = 1 : sample_size
                if scale == 1
                    if order_pre(t) <= ntrial(fid_seq(1))
                        trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_pre(t)+section_idx(fid_seq(1),1) - trial_len/2);
                    else
                        trial_flash_pre = binned(trial_len*(order_pre(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                            trial_len*(order_pre(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2);
                    end
                elseif scale == 2
                    trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_pre(t)+section_idx(fid_seq(1),1) - trial_len/2);
                end
                other_flash_pre = sum_flash_pre - trial_flash_pre;
                mean_flash_pre = other_flash_pre ./ (length(trial_num_pre) - 1) - mean_all .* length(trial_num_pre) ./ (length(trial_num_pre) - 1);

                if scale == 1
                    if order_post(t) <= ntrial(fid_seq(1))
                        trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 :...
                            trial_len*order_post(t)+section_idx(fid_seq(1),1));
                    else
                        trial_flash_post = binned(trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2 + 1 : ...
                            trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
                    end
                elseif scale == 2 
                    trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 : ...
                            trial_len*order_post(t)+section_idx(fid_seq(1),1));
                end
                other_flash_post = sum_flash_post - trial_flash_post;
                mean_flash_post = other_flash_post ./ (length(trial_num_post) - 1) - mean_all .* length(trial_num_post) ./ (length(trial_num_post) - 1);
                
                discriminant = (mean_flash_pre - mean_flash_post)';
                corrpos(t) = (trial_flash_pre - trial_flash_post) * discriminant; % no need to zero mean trial_flash & _null because they cancel out
            end
            corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
            Pc(flash_intensity - 1, test) = corr / length(corrpos);
        end
    end

    Pc_avg = mean(Pc,2);
    Pc_var = std(Pc,1,2);
    
%     if max(Pc_avg) >= 0.84
%         x = 1 : length(marker)-3; % exclude dark & 1ms flash (990, 52.1, 42.1)
%         flash_1ms = abs(marker(2:end) - floor(marker(2:end)))-0.1 < 1e-4;
        Pc_avg = Pc_avg(~flash_1ms);
        Pc_var = Pc_var(~flash_1ms);
        c = x;
        r = Pc_avg;
        fit = fitNakaRushton(c,r)

        hold on
        c = linspace(min(x), max(x), 10^3);
        Pc_fit = fit.Rmax * ((c.^fit.n) ./ ((c.^fit.n) + fit.c50.^fit.n)) + fit.offset;
        line([min(x), max(x)], [1, 1],'Color', [0 1 0])
        line([min(x), max(x)], [0.84, 0.84],'Color', [0 1 0])
        plot(c, Pc_fit, 'b--', 'LineWidth', 1)
%         scatter(c, Pc_avg, 'r')
        errorbar(x, Pc_avg, Pc_var, 'Color', [1 0 0])

        ylim([0.4, 1.05])
        titleStr = sprintf('c50: %0.2f n: %0.2f\n  Rmax: %0.2f offset: %0.2f', fit.c50,fit.n, fit.Rmax,fit.offset);
        title(sprintf(titleStr));
        xlabel('log(intensity)')
        ylabel('probability correct')
        ylim([0.4, 1.05])
        saveas(gcf, [flag num2str(ds_slave_now(cellnum)) '-2AFC-fit' '.png'])

%         x = 1 : length(marker)-1;
%         yline(1,'-.g'); yline(0.84,'-.g');
%         xticks(x)
%         xticklabels({'52.1','52.2','52.4','52.8','42.1','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})
%         xticklabels({'52.2','52.4','52.8','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})        
%         xtickangle(45)
                
%         saveas(gcf, [num2str(ds_slave_normal(cellnum)), '-2AFC-post_as_null', '.png'])
% %         print(['log_intensity-', num2str(ds_slave_normal(cellnum))], '-dpdf', '-fillpage')
        disp(['saved fig for ', num2str(ds_slave_now(cellnum))])
        close
%     else
%         cell_excluded = cell_excluded + 1;
%         disp([num2str(ds_slave_id_seq(cellnum)),' excluded due to low Pc'])
% %         close
%     end
end
disp('done')

%% fit single cell Naka Rushton
c = x;
r = Pc_avg;
fit = fitNakaRushton(c,r)

hold on
Pc_fit = fit.Rmax * ((c.^fit.n) ./ ((c.^fit.n) + fit.c50.^fit.n)) + fit.offset;
line([min(x), max(x)], [1, 1],'Color', [0 1 0])
line([min(x), max(x)], [0.84, 0.84],'Color', [0 1 0])
plot(c, Pc_fit, 'b--')
scatter(c, Pc_avg, 'r')

ylim([0.4, 1.05])
titleStr = sprintf('c50: %0.2f n: %0.2f\n  Rmax: %0.2f offset: %0.2f', fit.c50,fit.n, fit.Rmax,fit.offset);
title(sprintf(titleStr));
saveas(gcf, [num2str(ds_slave_normal(2)), '-2AFC-fit', '.png'])

%% solve detection threshold | BUG HERE
syms th
eqn = fit.Rmax * ((th.^fit.n) ./ ((th.^fit.n) + fit.c50.^fit.n)) + fit.offset - 0.84 == 0;
% S = simplify(solve(eqn))
% S = solve(eqn,th,'Real',true)
V = vpasolve(eqn,th,[min(c) max(c)])

%% flash vs null trial
% 
% % cell_excluded = 0;
% for c = length(ds_slave_normal)-1 : length(ds_slave_normal)
%     figure
%     ds_slave_index = find(datarun.cell_ids == ds_slave_normal(c)); 
%     spike_time = datarun.spikes{ds_slave_index, 1};
%     [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin
%     
%     sum_null = zeros(ntrial(1), trial_len); 
%     for t = 1 : ntrial(1)
%         trial_null = binned(trial_len*(t-1)+section_idx(1,1)+1 : trial_len*t+section_idx(1,1));
%         sum_null(t,:) = trial_null;
%     end
%     sum_null = sum(sum_null,1);
% 
%     ntest = 1000;
%     Pc = zeros(length(marker)-1 ,ntest);
%     for test = 1 : ntest
%         for flash_intensity = 2 : length(marker) % exclude dark==990. should improve by excluding 1ms here
%             nid = 1; 
%             fid = section_idx(:,end)==marker(flash_intensity);
%             fid_seq = find(fid==1);
%             
%             if floor(section_idx(fid_seq(1),7)) == 4
%                 scale = 2; % account for 4s trials
%             else
%                 scale = 1;
%             end
% 
%             sum_flash_seq = zeros(1, trial_len*scale);
%             for i = 1 : length(fid_seq)
%                 sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len*scale);
%                 for t = 1 : ntrial(fid_seq(i))
%                     trial_flash = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
%                     sum_flash_section{i}(t,:) = trial_flash;
%                 end
%                 sum_flash_section{i} = sum(sum_flash_section{i},1);
%                 sum_flash_seq = sum_flash_seq + sum_flash_section{i};
%             end
%             sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
%             
%             sum_all = sum_flash_seq + sum_null;
%             trial_num_null = 1 : ntrial(nid);
%             trial_num_flash = 1 : scale : scale*sum(ntrial(fid));
%             mean_all = sum_all ./ (length(trial_num_null) + length(trial_num_flash));
%             
%             sample_size = min(length(trial_num_null), length(trial_num_flash));
%             order_null = datasample(trial_num_null, sample_size, 'Replace', false);
%             order_flash = datasample(trial_num_flash, sample_size, 'Replace', false);
% 
%             corrpos = zeros(sample_size, 1);
%             for t = 1 : sample_size
%                 trial_null = binned(trial_len*(order_null(t)-1)+section_idx(nid,1)+1 : trial_len*order_null(t)+section_idx(nid,1));
%                 other_null = sum_null - trial_null;
%                 mean_null = other_null ./ (length(trial_num_null) - 1) - mean_all .* length(trial_num_null) ./ (length(trial_num_null) - 1); % zero-mean
% 
%                 if scale == 1
%                     if order_flash(t) <= ntrial(fid_seq(1))
%                         trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
%                             trial_len*order_flash(t)+section_idx(fid_seq(1),1));
%                     else
%                         trial_flash = binned(trial_len*(order_flash(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
%                             trial_len*(order_flash(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
%                     end
%                 elseif scale == 2
%                     trial_flash = binned(trial_len*(order_flash(t)-1)+section_idx(fid_seq(1),1)+1 : ...
%                             trial_len*order_flash(t)+section_idx(fid_seq(1),1));
%                 end
%                 other_flash = sum_flash_seq - trial_flash;
%                 mean_flash = other_flash ./ (length(trial_num_flash) - 1) - mean_all .* length(trial_num_flash) ./ (length(trial_num_flash) - 1);
% 
%                 discriminant = (mean_flash - mean_null)';
%                 corrpos(t) = (trial_flash - trial_null) * discriminant; % no need to zero mean trial_flash & _null because they cancel out
%             end
%             corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
%             Pc(flash_intensity - 1, test) = corr / length(corrpos);
%         end
%     end
% 
%     Pc_avg = mean(Pc,2);
%     Pc_var = std(Pc,1,2);
%     
% %     if max(Pc_avg) >= 0.84
% %         x = 1 : length(marker)-3; % exclude dark & 1ms flash (990, 52.1, 42.1)
% %         flash_1ms = abs(marker(2:end) - floor(marker(2:end)))-0.1 < 1e-4;
%         Pc_avg = Pc_avg(~flash_1ms);
%         Pc_var = Pc_var(~flash_1ms);
% 
% %         x = 1 : length(marker)-1;
%         errorbar(x, Pc_avg, Pc_var)
%         hold on
%         line([min(x), max(x)], [1, 1],'Color', [0 1 0])
%         line([min(x), max(x)], [0.84, 0.84],'Color', [0 1 0])
% %         yline(1,'-.g'); yline(0.84,'-.g');
% %         xticks(x)
% %         xticklabels({'52.1','52.2','52.4','52.8','42.1','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})
% %         xticklabels({'52.2','52.4','52.8','42.2','42.4','42.8','32.2','32.4','34.8','24.2'})        
% %         xtickangle(45)
%         xlabel('log(intensity)')
%         ylabel('probability correct')
%         ylim([0.4, 1.05])
%         
%         saveas(gcf, [num2str(ds_slave_normal(c)), '-2AFC.png'])
% %         print(['log_intensity-', num2str(ds_slave_normal(c))], '-dpdf', '-fillpage')
%         disp(['saved fig for ', num2str(ds_slave_normal(c))])
%         close
% %     else
% %         cell_excluded = cell_excluded + 1;
% %         disp([num2str(ds_slave_normal(c)),' excluded due to low Pc'])
% %         close
% %     end
% end
% disp('done')


%% test Discriminant shape for low performance cells

% low_perf = sort(ds_slave_normal);
% low_perf = low_perf(5:end);

for c = 1 : length(ds_slave_normal)

    ds_slave_index = find(datarun.cell_ids == ds_slave_normal(c)); 
    spike_time = datarun.spikes{ds_slave_index, 1};
    [binned, ~] = histcounts(spike_time, edges); % binned = vector of nspike in each 20 ms bin
    D_all = cell(length(marker),1);

    ntest = 1000;
    Pc = zeros(length(marker)-1 ,ntest);
    for test = 1 : ntest
        for flash_intensity = 2 : length(marker) % exclude dark==990. should improve by excluding 1ms here
%             nid = 1; 
            fid = section_idx(:,end)==marker(flash_intensity);
            fid_seq = find(fid==1);
            D_fid = zeros(ntest*120, 50);
            
            if floor(section_idx(fid_seq(1),7)) == 4
                scale = 2; % account for 4s trials
            else
                scale = 1;
            end

            sum_flash_seq = zeros(1, trial_len*scale);
            for i = 1 : length(fid_seq)
                sum_flash_section{i} = zeros(ntrial(fid_seq(i)), trial_len*scale);
                for t = 1 : ntrial(fid_seq(i))
                    trial_flash_full = binned(trial_len*scale*(t-1)+section_idx(fid_seq(i),1)+1 : trial_len*scale*t+section_idx(fid_seq(i),1));
                    sum_flash_section{i}(t,:) = trial_flash_full;
                end
                sum_flash_section{i} = sum(sum_flash_section{i},1);
                sum_flash_seq = sum_flash_seq + sum_flash_section{i};
            end
            sum_flash_seq = sum_flash_seq(1 : trial_len); % take only 0-2s of 4s trials
            sum_flash_pre = sum_flash_seq(1 : length(sum_flash_seq)/2); % pre (0-1s) as flash trial
            sum_flash_post = sum_flash_seq(length(sum_flash_seq)/2 + 1 : end); % post (1-2s) as null trial
            
            sum_all = sum_flash_pre + sum_flash_post;
            trial_num_post = 1 : scale : scale*sum(ntrial(fid));
            trial_num_pre = 1 : scale : scale*sum(ntrial(fid));
            mean_all = sum_all ./ (length(trial_num_post) + length(trial_num_pre)); 
            
            sample_size = min(length(trial_num_post), length(trial_num_pre));
            order_post = datasample(trial_num_post, sample_size, 'Replace', false);
            order_pre = datasample(trial_num_pre, sample_size, 'Replace', false);

            corrpos = zeros(sample_size, 1);
            for t = 1 : sample_size
                if scale == 1
                    if order_pre(t) <= ntrial(fid_seq(1))
                        trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_pre(t)+section_idx(fid_seq(1),1) - trial_len/2);
                    else
                        trial_flash_pre = binned(trial_len*(order_pre(t)-ntrial(fid_seq(1))-1)+section_idx(fid_seq(2),1)+1 : ...
                            trial_len*(order_pre(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2);
                    end
                elseif scale == 2
                    trial_flash_pre = binned(trial_len*(order_pre(t)-1)+section_idx(fid_seq(1),1)+1 : ...
                            trial_len*order_pre(t)+section_idx(fid_seq(1),1) - trial_len/2);
                end
                other_flash_pre = sum_flash_pre - trial_flash_pre;
                mean_flash_pre = other_flash_pre ./ (length(trial_num_pre) - 1) - mean_all .* length(trial_num_pre) ./ (length(trial_num_pre) - 1);

                if scale == 1
                    if order_post(t) <= ntrial(fid_seq(1))
                        trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 :...
                            trial_len*order_post(t)+section_idx(fid_seq(1),1));
                    else
                        trial_flash_post = binned(trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2 + 1 : ...
                            trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
                    end
                elseif scale == 2 
                    trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 : ...
                            trial_len*order_post(t)+section_idx(fid_seq(1),1));
                end
                other_flash_post = sum_flash_post - trial_flash_post;
                mean_flash_post = other_flash_post ./ (length(trial_num_post) - 1) - mean_all .* length(trial_num_post) ./ (length(trial_num_post) - 1);
                
                discriminant = (mean_flash_pre - mean_flash_post)';
                D_fid((test-1)*60 + t, :) = discriminant;
                corrpos(t) = (trial_flash_pre - trial_flash_post) * discriminant; % no need to zero mean trial_flash & _null because they cancel out
            end
            corr = sum(corrpos>0) + 1/2 * sum(corrpos==0);
            Pc(flash_intensity - 1, test) = corr / length(corrpos);
            D_all{flash_intensity,1} = mean(D_fid,1);

        end
    end

    Pc_avg = mean(Pc,2);
    Pc_var = std(Pc,1,2);
    Pc_avg = Pc_avg(~flash_1ms);
    Pc_var = Pc_var(~flash_1ms);
    
    intensity_seq = [3,4,5,7,8,9,10,11,12,13]; % exclude dark & 1ms flash
    % intensity_seq = 1:10;
    figure('units','normalized','outerposition',[0 0 1 1]) 
    for i = 1:10
        subplot(10,1,i)
        D_intensity = D_all{intensity_seq(i),1};
        plot(D_intensity)
        set(gca,'XTick',[], 'YTick', [])
    end
    saveas(gcf, [num2str(ds_slave_normal(c)), '-discriminant-exclude1ms-postasnull.png'])
    close
    disp(['saved fig for ', num2str(ds_slave_normal(c))])
end
disp('done')
