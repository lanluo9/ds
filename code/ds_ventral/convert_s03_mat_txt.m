
clear 
clc

load('2016-03-04_s03.mat')

%%
heading_str = '(:class :MG :spatial_modulation :square :frames 480 :x_start 0 :x_end 800 :y_start 0 :y_end 600)\n';
line_template = '(:DIRECTION 225 :TEMPORAL_PERIOD 120 :SPATIAL_PERIOD 240 :RGB #(0.48 0.48 0.48) :BACK_RGB #(0.5 0.5 0.5))';

%% 

fileID = fopen('s03.txt','w');
fprintf(fileID, heading_str);

line_end_str = ' :RGB #(0.50 0.50 0.50) :BACK_RGB #(0.25 0.25 0.25))\n';
sp = stim_out.spatial_period;

for i = 1 : length(stim_out.trials)
    dir = stim_out.trials(i).direction;
    tp = 10 * stim_out.trials(i).temporal_period; % tp naming convention btw 2016 vs 2019? 2019 = 120, 240
    fprintf(fileID, ['(:DIRECTION ', num2str(dir), ' :TEMPORAL_PERIOD ', num2str(tp), ' :SPATIAL_PERIOD ', num2str(sp), line_end_str]);
    
end

fclose(fileID);
