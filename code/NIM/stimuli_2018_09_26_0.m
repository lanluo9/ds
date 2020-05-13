%% 2018-09-26-0
% JC

%% Initialization

my_path = '/Users/acquisition/Duke-devo/';

addpath(genpath(my_path))
cd(my_path)

path2save = [my_path, '/saved_timestamps/2016-08-03/'];
screen_number = 2; %0=test on current monitor, 2= is two monitor display
def_params = initialize_display('OLED', screen_number); % default parameters
%mglMoveWindow([])

% real refresh rate 
%mglTestRefresh(2)

% set gamma OLED Aug 2 2016
scale = [1.1399    1.0998    1.1027];
power = [1.1741    1.2998    1.3112];
offset = [-0.1445   -0.1023   -0.1054];
set_gamma_from_fit_params(scale, power, offset);


%% Focus Squares

fprintf('\n\n<strong> Focus squares. </strong>\n');
clear parameters stimulus;

parameters.class = 'FS';
stimulus = make_stimulus(parameters,def_params);
display_stimulus(stimulus);


%% black
mglClearScreen(0);
mglFlush

%% DS typing with full tuning curves

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
                                           'TEMPORAL_PERIOD', [12 24 60 120 240 480 720 960 1440], ...
                                           'SPATIAL_PERIOD', [240], ...
                                           'RGB', {[0.25 0.25 0.25]}, ...
                                           'BACK_RGB', {[0.5 0.5 0.5]}, ...
                                           'nrepeats',4);
path2file = write_s_file(parameters, variable_parameters);
s_params = read_s_file(path2file);

% see second option example in "S File read"
for i=2:size(s_params,2)
    trial_params = combine_parameters(s_params{1}, s_params{i});
    stimulus{i-1} = make_stimulus(trial_params, def_params);
end
for i=1:length(stimulus)
    display_stimulus(stimulus{i}, 'wait_trigger', 0, 'wait_key', 0);
end

%% Image Jitter - repeat

parameters.class = 'IJ'; 
 
parameters.image_path = '/Users/acquisition/Desktop/Natural_Images/VhNaturalImages/VhImage_imk03002.mat' ; % (path of background image to jitter)
parameters.image_display_flag = true ; % (true = display background image)
parameters.square_display_flag = false ; % (true = display square ontop of image)               
parameters.x_start = 100 ; % mask bounds
parameters.x_end = 500 ;
parameters.y_start = 200 ;
parameters.y_end =  600 ;   
 
parameters.image_jitter_std = 60 ; % (pix) std of image x,y shifts
parameters.square_jitter_std = 0 ; % (pix) std of square x,y shifts
parameters.jitter_smooth_frames = 300 ; % (frame number) over which position is averaged
parameters.square_width = 100 ; % width of square
parameters.square_intensity_std = -6 ; % (std) of image of the square intensity (- is a decrement)
parameters.num_repeats = 20 ;
parameters.jitter_seed = 11111 ;
parameters.num_frames = 60*10 ; % number of independent stimuli frames (time = (num_frames*frame_interval)/frame_rate)
parameters.frame_interval = 1 ; % number of frames each stim frame is displayed (refresh)
                
parameters.wait_trigger = 0;                            
parameters.wait_key = 0 ; 
 
stimulus = make_stimulus(def_params, parameters);
display_stimulus(stimulus,'wait_trigger',parameters.wait_trigger, 'wait_key',parameters.wait_key) ;

%% mounted

%% black
mglClearScreen(0.5);
mglFlush

%% BW focus
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
parameters.x_start = 101;  parameters.x_end = 700;
parameters.y_start = 1;   parameters.y_end = 600;

parameters.independent = 0;
parameters.interval = 6;
parameters.stixel_width = 5;
parameters.frames = 60000;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 0,'erase',0);

%% data000 DS typing with full speed tuning curves

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
                                           'TEMPORAL_PERIOD', [12 24 60 120 240 480 720 960 1440], ...
                                           'SPATIAL_PERIOD', [240], ...
                                           'RGB', {[0.25 0.25 0.25]}, ...
                                           'BACK_RGB', {[0.5 0.5 0.5]}, ...
                                           'nrepeats',4);
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


%% data001 Image Jitter

parameters.class = 'IJ'; 
 
parameters.image_path = '/Users/acquisition/Desktop/Natural_Images/VhNaturalImages/VhImage_imk03002.mat' ; % (path of background image to jitter)
parameters.image_display_flag = true ; % (true = display background image)
parameters.square_display_flag = false ; % (true = display square ontop of image)               
parameters.x_start = 100 ; % mask bounds
parameters.x_end = 500 ;
parameters.y_start = 200 ;
parameters.y_end =  600 ;   
 
parameters.image_jitter_std = 60 ; % (pix) std of image x,y shifts
parameters.square_jitter_std = 0 ; % (pix) std of square x,y shifts
parameters.jitter_smooth_frames = 300 ; % (frame number) over which position is averaged
parameters.square_width = 100 ; % width of square
parameters.square_intensity_std = -6 ; % (std) of image of the square intensity (- is a decrement)
parameters.num_repeats = 1 ;
parameters.jitter_seed = 11111 ;
parameters.num_frames = 60*60*30 ; % number of independent stimuli frames (time = (num_frames*frame_interval)/frame_rate)
parameters.frame_interval = 1 ; % number of frames each stim frame is displayed (refresh)
                
parameters.wait_trigger = 1;                            
parameters.wait_key = 0 ; 
 
stimulus = make_stimulus(def_params, parameters);
display_stimulus(stimulus,'wait_trigger',parameters.wait_trigger, 'wait_key',parameters.wait_key) ;


%% datat002 ds typing - short

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


%% data003 Image Jitter

parameters.class = 'IJ'; 
 
parameters.image_path = '/Users/acquisition/Desktop/Natural_Images/VhNaturalImages/VhImage_imk00001.mat' ; % (path of background image to jitter)
parameters.image_display_flag = true ; % (true = display background image)
parameters.square_display_flag = false ; % (true = display square ontop of image)               
parameters.x_start = 100 ; % mask bounds
parameters.x_end = 500 ;
parameters.y_start = 200 ;
parameters.y_end =  600 ;   
 
parameters.image_jitter_std = 60 ; % (pix) std of image x,y shifts
parameters.square_jitter_std = 0 ; % (pix) std of square x,y shifts
parameters.jitter_smooth_frames = 300 ; % (frame number) over which position is averaged
parameters.square_width = 100 ; % width of square
parameters.square_intensity_std = -6 ; % (std) of image of the square intensity (- is a decrement)
parameters.num_repeats = 1 ;
parameters.jitter_seed = 11111 ;
parameters.num_frames = 60*60*30 ; % number of independent stimuli frames (time = (num_frames*frame_interval)/frame_rate)
parameters.frame_interval = 1 ; % number of frames each stim frame is displayed (refresh)
                
parameters.wait_trigger = 1;                            
parameters.wait_key = 0 ; 
 
stimulus = make_stimulus(def_params, parameters);
display_stimulus(stimulus,'wait_trigger',parameters.wait_trigger, 'wait_key',parameters.wait_key) ;


%% data004 Image Jitter - repeat

parameters.class = 'IJ'; 
 
parameters.image_path = '/Users/acquisition/Desktop/Natural_Images/VhNaturalImages/VhImage_imk03002.mat' ; % (path of background image to jitter)
parameters.image_display_flag = true ; % (true = display background image)
parameters.square_display_flag = false ; % (true = display square ontop of image)               
parameters.x_start = 100 ; % mask bounds
parameters.x_end = 500 ;
parameters.y_start = 200 ;
parameters.y_end =  600 ;   
 
parameters.image_jitter_std = 60 ; % (pix) std of image x,y shifts
parameters.square_jitter_std = 0 ; % (pix) std of square x,y shifts
parameters.jitter_smooth_frames = 300 ; % (frame number) over which position is averaged
parameters.square_width = 100 ; % width of square
parameters.square_intensity_std = -6 ; % (std) of image of the square intensity (- is a decrement)
parameters.num_repeats = 20 ;
parameters.jitter_seed = 11111 ;
parameters.num_frames = 40*100 ; % number of independent stimuli frames (time = (num_frames*frame_interval)/frame_rate)
parameters.frame_interval = 1 ; % number of frames each stim frame is displayed (refresh)
                
parameters.wait_trigger = 1 ;                            
parameters.wait_key = 0 ; 
 
stimulus = make_stimulus(def_params, parameters);
display_stimulus(stimulus,'wait_trigger',parameters.wait_trigger, 'wait_key',parameters.wait_key) ;



%% data005 
% flow without jitter

%% data006
% flow with jitter

%% data007
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
% parameters.x_start = 1;  parameters.x_end = 800;
% parameters.y_start = 1;   parameters.y_end = 600;

%%%%%%%%%%%%%% CRT %%%%%%%%%%%%%% 
parameters.x_start = 2;  parameters.x_end = 797;

parameters.y_start = 0;   parameters.y_end = 600 ;

parameters.independent = 0;
parameters.interval = 1;
parameters.stixel_width = 15;
parameters.frames = 60*60*60;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

% For Voronoi, set stixel_height and stixel_width to 1 and pass a map path
%parameters.map_file_name = [my_path, '/Maps/2011-12-13-2_f04_vorcones/map-0000.txt'];
% parameters.map_file_name = ['/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Stimulus Code/test/data002/large_on/5.txt'];

% parameters.map_file_name = ['/Volumes/Data/2016-01-05-1/Visual/maps/map_data001.txt'];
% mask example
% mask = zeros(80,40);
% mask(1:10, 1:10) = 255;
% parameters.mask = mask;

% parameters.map_file_name = '/Volumes/Lab/Users/crhoades/Colleen/matlab/private/colleen/New Cell Types/Stimulus Code/test/data002/large_on/5.txt';
stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 1, 'erase', 1);


%% BW focus
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
parameters.x_start = 101;  parameters.x_end = 700;
parameters.y_start = 1;   parameters.y_end = 600;

parameters.independent = 0;
parameters.interval = 6;
parameters.stixel_width = 50;
parameters.frames = 1;

parameters.stixel_height = parameters.stixel_width;
parameters.field_width = (parameters.x_end-parameters.x_start+1)/parameters.stixel_width;  
parameters.field_height = (parameters.y_end-parameters.y_start+1)/parameters.stixel_height;

stimulus = make_stimulus(parameters, def_params);

time_stamps = display_stimulus(stimulus, 'wait_trigger', 0,'erase',0);


%% end photons
Stop_Photons

make_normal_table(8)
