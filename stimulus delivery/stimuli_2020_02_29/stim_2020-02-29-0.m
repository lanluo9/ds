%%Greg 2020-02-29-0 C57-BL6

%% Initialization

my_path = '/Users/stimulus/Desktop/Duke-devo';

addpath(genpath(my_path))
cd(my_path)

path2save = [my_path, '/saved_stim/2016-08-03/'];
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

%% data000
% flashes from the function generator and LED


%% data001 Moving Grating S File write
%NDF3
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
                                           'TEMPORAL_PERIOD', [120 240 480], ...
                                           'SPATIAL_PERIOD', [240], ...
                                           'RGB', {[0.48 0.48 0.48]}, ...
                                           'BACK_RGB', {[0.5 0.5 0.5]}, ...
                                           'nrepeats',6);
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

%% data002 Moving Grating S File write
%NDF 2
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
                                           'TEMPORAL_PERIOD', [120 240 480], ...
                                           'SPATIAL_PERIOD', [240], ...
                                           'RGB', {[0.48 0.48 0.48]}, ...
                                           'BACK_RGB', {[0.5 0.5 0.5]}, ...
                                           'nrepeats',6);
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

%% data003 Random Noise

fprintf('\n\n<strong> Random Noise </strong>\n');
clear parameters stimulus

parameters.class = 'RN';
parameters.back_rgb = [1 1 1]*0.5;
parameters.rgb = [1 1 1]*0.48;
parameters.seed = 11111;
parameters.binary = 1;
parameters.probability = 1;
parameters.jitter = 0;
parameters.delay_frames = 0;

%%%%%%%%%%%%%% OLED %%%%%%%%%%%%%% 
parameters.x_start = 1;  parameters.x_end = 800;
parameters.y_start = 1;   parameters.y_end = 600;

%%%%%%%%%%%%%% CRT %%%%%%%%%%%%%% 
% parameters.x_start = 1;  parameters.x_end = 640;
% parameters.y_start = 1;   parameters.y_end = 480;

parameters.independent = 0;
parameters.interval = 3;
parameters.stixel_width = 10;
parameters.frames = 60*3600;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height

% For Voronoi, set stixel_height and stixel_width to 1 and pass a map path
% parameters.map_file_name = [my_path, '/Maps/2011-12-13-2_f04_vorcones/map-0000.txt'];

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 1);




