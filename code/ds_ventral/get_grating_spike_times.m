                  function spike_times = get_grating_spike_times(datarun, cell_spec, stimulus_diration, varargin)
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
% direction       first dir         gets first direction in datarun.stimulus
% TP              first TP          gets first TP in datarun.stimulus
% SP              first SP         	gets first SP in datarun.stimulus
% contrast        []                only extracts if specified (may cause error if not specified)
% background      []                only extracts if specified
%
%
% Created: GDF, 2019-08-15
%

% BEGIN FUNCTION
% sort varargin
p = inputParser;
p.addParameter('direction', datarun.stimulus.params.DIRECTION(1), @isnumeric);
p.addParameter('TP', datarun.stimulus.params.TEMPORAL_PERIOD(1), @isnumeric);
p.addParameter('SP', datarun.stimulus.params.SPATIAL_PERIOD(1), @isnumeric);
p.addParameter('contrast', []);
p.addParameter('background', []);

p.parse(varargin{:});

% get indices to specified RGC list or type
cell_indices = get_cell_indices(datarun, cell_spec);
num_rgcs = length(cell_indices);

% extract the number of times each stimulus was presented.
num_repeats = datarun.stimulus.repetitions;
num_stim = length(datarun.stimulus.combinations);


% initialize output cell of spike times
spike_times = cell(num_rgcs, num_repeats);

% find the numeric index to the stim matching the user's input
for gstim = 1:num_stim
    if datarun.stimulus.combinations(gstim).TEMPORAL_PERIOD == p.Results.TP && ...
        datarun.stimulus.combinations(gstim).DIRECTION == p.Results.direction && ...
        datarun.stimulus.combinations(gstim).SPATIAL_PERIOD == p.Results.SP
    stim_index = gstim;
    break
    end
end

% loop over rgcs and stim repeats to extract cell array of spike times
for rgc = 1:num_rgcs
    temp_trigger_inds = find(datarun.stimulus.trial_list == stim_index);
    
    if length(temp_trigger_inds) ~= num_repeats
        error('insufficient numbr of trials found to match the specified number of repeats');
    end
    
    for g_rep = 1:num_repeats
        tmp_spike_times = datarun.spikes{rgc}(datarun.spikes{rgc} >= datarun.stimulus.triggers(temp_trigger_inds(g_rep)) &...
                            datarun.spikes{rgc} < datarun.stimulus.triggers(temp_trigger_inds(g_rep)) + stimulus_diration);
        tmp_spike_times = tmp_spike_times - datarun.stimulus.triggers(temp_trigger_inds(g_rep));
        
        spike_times{rgc, g_rep} = tmp_spike_times;
    end
end



