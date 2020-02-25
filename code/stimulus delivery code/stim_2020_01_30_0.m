%%Miranda 2020-01-30 Cngb1-/- 3M
%% Initialization

my_path = '/Users/jcafaro/Documents/MATLAB/Photons';
%my_path = '/Users/acquisition/Photons';
my_path = '/Users/acquisition/Duke-devo/';

addpath(genpath(my_path))
cd(my_path)

path2save = [my_path, '/saved_timestamps/2016-08-03/'];
screen_number = 2; %0=test on current monitor, 2= is two monitor display
def_params = initialize_display('OLED', screen_number); % default parameters
%mglMoveWindow([])

% real refresh rate 
%mglTestRefresh(2)

% set gamma OLED Nov 14 2019
scale = [1.03568066871491,1.03647726974954,1.03801579515696];
power = [1.07164946067619,1.05540770812658,1.03021347204221];
offset = [-0.00771159700138625,-0.0207938158112673,-0.0249832730848717];
set_gamma_from_fit_params(scale, power, offset);

%% Focus Squares

fprintf('\n\n<strong> Focus squares. </strong>\n');
clear parameters stimulus;
 
parameters.class = 'FS';
stimulus = make_stimulus(parameters,def_params);
display_stimulus(stimulus);

%% Set background

% white
mglClearScreen(0);
mglFlush

%%data000 spont act in dark

%% Set background

% black
mglClearScreen(0.5);
mglFlush


%% data001 is just gray screen (spontaneous activity)

%%  white noise repeats data002
%%currently set to 10s BW with 200 repeats, can lower to 5s if short on
%%time

fprintf('\n\n<strong> Random Noise Repeats </strong>\n');
clear parameters stimulus time_stamps

num_rep=200;
rep_length = 10;
parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.48;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

%%%%%%%%%%%%%% OLED %%%%%%%%%%%%%% 
parameters.x_start = 10;  parameters.x_end = 790;
parameters.y_start = 0;   parameters.y_end = 600;

%%%%%%%%%%%%%% CRT %%%%%%%%%%%%%% 
% parameters.x_start = 1;  parameters.x_end = 640;
% parameters.y_start = 1;   parameters.y_end = 480;

parameters.independent = 0;
parameters.interval =4;
parameters.stixel_width = 30;
parameters.frames = (60*rep_length)-4;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);
for i=1:num_rep
    time_stamps{i} = display_stimulus(stimulus, 'wait_trigger',1, 'erase', 0);
end

%% Random Noise data003

fprintf('\n\n<strong> Random Noise </strong>\n');
clear parameters stimulus

parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.49;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

parameters.x_start = 10;  parameters.x_end = 790;

parameters.y_start = 0;   parameters.y_end = 600 ;

parameters.independent = 0;
parameters.interval = 4;
parameters.stixel_width = 30;
parameters.frames = 60*60*30;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 1, 'erase', 0);

%% Random Noise to focus

fprintf('\n\n<strong> Random Noise </strong>\n');
clear parameters stimulus

parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.49;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

parameters.x_start = 0;  parameters.x_end = 800;

parameters.y_start = 0;   parameters.y_end = 600 ;

parameters.independent = 0;
parameters.interval = 4;
parameters.stixel_width = 30;
parameters.frames = 60*60*30;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 0, 'erase', 0);

%%  white noise repeats data004
%%currently set to 10s BW with 200 repeats, can lower to 5s if short on
%%time

fprintf('\n\n<strong> Random Noise Repeats </strong>\n');
clear parameters stimulus time_stamps

num_rep=200;
rep_length = 10;
parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.48;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

%%%%%%%%%%%%%% OLED %%%%%%%%%%%%%% 
parameters.x_start = 10;  parameters.x_end = 790;
parameters.y_start = 0;   parameters.y_end = 600;

%%%%%%%%%%%%%% CRT %%%%%%%%%%%%%% 
% parameters.x_start = 1;  parameters.x_end = 640;
% parameters.y_start = 1;   parameters.y_end = 480;

parameters.independent = 0;
parameters.interval =2;
parameters.stixel_width = 15;
parameters.frames = (60*rep_length)-4;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);
for i=1:num_rep
    time_stamps{i} = display_stimulus(stimulus, 'wait_trigger',1, 'erase', 0);
end

%% Random Noise data005

fprintf('\n\n<strong> Random Noise </strong>\n');
clear parameters stimulus

parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.49;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

parameters.x_start = 10;  parameters.x_end = 790;

parameters.y_start = 0;   parameters.y_end = 600 ;

parameters.independent = 0;
parameters.interval = 2;
parameters.stixel_width = 15;
parameters.frames = 60*60*30;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 1, 'erase', 0);

%% data006 ds typing

fprintf('\n\n<strong> Moving Grating. </strong>\n');
clear parameters stimulus

parameters.class = 'MG';
parameters.spatial_modulation = 'square'; % sine or square
% parameters.rgb = [1 1 1]*0.25;
% parameters.back_rgb = [1 1 1]*0.5;
parameters.frames = 8*60; % presentation of each grating, frames
parameters.x_start = 0;  parameters.x_end = 800;
parameters.y_start = 0;   parameters.y_end = 600;
% parameters.direction = 45;

variable_parameters = randomize_parameters('DIRECTION', [0 45 90 135 180 225 270 315], ...
                                           'TEMPORAL_PERIOD', [60 240], ...
                                           'SPATIAL_PERIOD', [240], ...
                                           'RGB', {[0.49 0.49 0.49]}, ...
                                           'BACK_RGB', {[0.5 0.5 0.5]}, ...
                                           'nrepeats',5);
path2file = write_s_file(parameters, variable_parameters);
s_params = read_s_file(path2file);

% see second option example in "S File read"
for i=2:size(s_params,2)
    trial_params = combine_parameters(s_params{1}, s_params{i});
    stimulus{i-1} = make_stimulus(trial_params, def_params);
end
for i=1:length(stimulus)
    display_stimulus(stimulus{i}, 'wait_trigger', 1, 'wait_key', 0);
end

%% natural repeats - cat data007

fprintf('\n\n<strong> Raw Movie </strong>\n');
clear parameters stimulus;

num_rep=90;
rep_length = 10;
parameters.class = 'RM';
parameters.back_rgb = [1 1 1]*0.5;
parameters.x_start = 80; %80 if change oled % x_end and y_end wil depend on movie size (and stixel size)!
parameters.y_start = 60;
parameters.stixel_width = 2;   parameters.stixel_height = 2;
parameters.frames = (60*rep_length)-4;  % duration of each repetition, default whole movie
parameters.start_frame = 0; % >0
parameters.interval = 1;
parameters.flip = 1;  % 1 = normal; 2 = vertical flip; 3 = horizontal flip; 4 = vertical + horizontal flip
parameters.reverse = 0;   % 1 = backward (reverse), 0 = forward

parameters.movie_name = '/Users/acquisition/Desktop/Movies/cat_mean117_sd62_0to255.rawMovie';

stimulus = make_stimulus(parameters, def_params);

for i = 1:num_rep
    time_stamps = display_stimulus(stimulus,'trigger_interval', 100,'wait_key',0, 'erase', 0,'wait_trigger', 1);
end

%% natural repeats - squirrel data008

fprintf('\n\n<strong> Raw Movie </strong>\n');
clear parameters stimulus;

num_rep=90;
rep_length = 10;
parameters.class = 'RM';
parameters.back_rgb = [1 1 1]*0.5;
parameters.x_start = 80; %80 if change oled % x_end and y_end wil depend on movie size (and stixel size)!
parameters.y_start = 60;
parameters.stixel_width = 2;   parameters.stixel_height = 2;
parameters.frames = (60*rep_length)-4;  % duration of each repetition, default whole movie
parameters.start_frame = 0; % >0
parameters.interval = 4;
parameters.flip = 1;  % 1 = normal; 2 = vertical flip; 3 = horizontal flip; 4 = vertical + horizontal flip
parameters.reverse = 0;   % 1 = backward (reverse), 0 = forward

parameters.movie_name = '/Users/acquisition/Desktop/Movies/squirrel_mean117_sd62_0to255.rawMovie';

stimulus = make_stimulus(parameters, def_params);

for i = 1:num_rep
    time_stamps = display_stimulus(stimulus,'trigger_interval', 100,'wait_key',0, 'erase', 0,'wait_trigger', 1);
end