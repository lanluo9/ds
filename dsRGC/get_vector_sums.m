function [vector_sums, vector_mags] = get_vector_sums(datarun, cell_spec, varargin)
%
% usage: function spike_times = get_grating_spike_times(datarun, cell_ids, varargin)
%
% arguments:  datarun - datarun struct
%           cell_spec - which cells (see get_cell_indices for options)
%           stimulus_duration - duration of the stimulus (e.g. typically 8 s).
%            varargin - struct or list of optional parameters (see below)
%
% outputs:    spike_times - a cell array of spike times that is num_rgcs x
%                           num_repeats
%
% optional parameters, their default values, and what they specify:
%
%
% TP              first TP          gets first TP in datarun.stimulus
% SP              first SP         	gets first SP in datarun.stimulus
% stim_duration          8               duration of the grating
% contrast        []                only extracts if specified (may cause error if not specified)
% background      []                only extracts if specified
%
%
% Created: GDF, 2019-08-15
%

p = inputParser;
p.addParameter('TP', datarun.stimulus.params.TEMPORAL_PERIOD(1), @isnumeric);
p.addParameter('SP', datarun.stimulus.params.SPATIAL_PERIOD(1), @isnumeric);
p.addParameter('stim_duration', 8, @isnumeric);
p.addParameter('contrast', []);
p.addParameter('background', []);
p.parse(varargin{:});

cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);

%vector_sums = zeros(num_rgcs, 2);
%vector_mags = zeros(num_rgcs, 1);
spike_counts = zeros(num_rgcs, length(datarun.stimulus.params.DIRECTION));

for g_dir = 1:length(datarun.stimulus.params.DIRECTION)
        temp_dir = datarun.stimulus.params.DIRECTION(g_dir);
        tmp_spike_times = get_grating_spike_times(datarun, cell_spec, p.Results.stim_duration, 'direction', temp_dir, 'SP', p.Results.SP, 'TP', p.Results.TP);
        tmp_gratingrun.direction(g_dir).spike_times = tmp_spike_times;
end

for rgc = 1:num_rgcs
    for g_dir = 1:length(datarun.stimulus.params.DIRECTION)
        tmp_counter = 0;
        for g_rep = 1:datarun.stimulus.repetitions
            tmp_counter = tmp_counter + length(tmp_gratingrun.direction(g_dir).spike_times{rgc, g_rep});
        end
        spike_counts(rgc, g_dir) = tmp_counter;
    end
end


for rgc = 1:num_rgcs
    norm_factor = max(spike_counts(rgc, :));
    
    x_cords = (spike_counts(rgc,:) ./ norm_factor) .* sind(datarun.stimulus.params.DIRECTION);
    y_cords = (spike_counts(rgc,:) ./ norm_factor) .* cosd(datarun.stimulus.params.DIRECTION);

    vector_sums(rgc,:) = [sum(x_cords), sum(y_cords)];
    vector_mags(rgc) = norm([sum(x_cords), sum(y_cords)]);   
end

            


       
