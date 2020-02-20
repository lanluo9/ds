%% import data

time_white = importdata('whitenoise600k500us5hz3_5mVSD.xlsx'); 
time_seq = time_white(:,1);
white = time_white(:,2);

event = importdata('Events20206008.txt');
event = event(~isnan(event)); % delete NaNs

%% sanity checks

white_mean = mean(white)
white_var = std(white) 
% figure; plot(white)

if min(event) < min(time_seq) || max(event) > max(time_seq)
    fprintf('event time range does not match white noise time range\n')
end
% figure; scatter(1:length(event), event)

%% event triggering segment of white noise

ets_len = 1000; % define event triggering segment length as 1000 ms
ets_start = event - ets_len;
ets_start_end = [ets_start, event];

ets = zeros(length(ets_start_end), ets_len/0.5);
for i = 1:length(ets_start_end)
    ets(i, :) = white(white(:,1)==ets_start_end());
end