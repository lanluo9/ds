%% import data

clc
clear

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

time_step = 0.5; % voltage and white noise is recorded every 0.5 ms
ets_len = 1000; % define event triggering segment length as 1000 ms
ets_start_end = [event - ets_len, event];

ets = zeros(length(ets_start_end), ets_len/time_step);
for i = 1:length(ets_start_end)
    id1 = find(time_seq==ets_start_end(i,1)) + 1;
    id2 = find(time_seq==ets_start_end(i,2));
    ets(i, :) = white(id1:id2)'; % ets is a matrix whose every row is a white noise voltage seq preceding an event
end

%% event triggering average

eta = mean(ets,1);
time_axis = [(-1*ets_len+time_step) : time_step : 0];

plot(time_axis, eta, 'LineWidth',2)
hold on 
line([-1000,0],[0 0],'Color','Black')

title('Event triggering average')
xlabel('Time (ms)')
ylabel('Voltage (mV)')


