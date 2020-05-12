function [mov,height,width,duration,refresh] = get_movie_LL(mdf_file, triggers, frames)
% GET_MOVIE      Load a white noise movie
%
%   length(triggers) must be >= 10*stimulus_interval(1,2...)
%
%greschner
% 
% make changes to take less memory when initializing & take out the color dimension (binary white only)
% input: path_to_xml, triggers, nframes
% LL

% load movie
[mvi] = load_movie(mdf_file, triggers);

% grab movie parameters: height, width, duration
height = double(mvi.getHeight);
width = double(mvi.getWidth);
duration = double(mvi.size);
refresh = double(mvi.getRefreshTime);

if exist('frames','var')
    if frames<=duration
        tduration=frames;
    else
        error('movie too short');
    end
end
  

if length(tduration) > 1
    mov=zeros(height,width,1,length(tduration));
    cntr = 0;
    for i=tduration(1):tduration(length(tduration))
        F = mvi.getFrame(i-1).getBuffer;
        
        F = reshape(F,3,width,height);
%         sum(sum(F(1,:,:) == F(2,:,:))) == height*width % sanity check: see if movie is really grayscale
%         sum(sum(F(1,:,:) == F(3,:,:))) == height*width
        F = F(1,:,:);
        cntr = cntr + 1;
        mov(:,:,:,cntr) = permute(F,[3 2 1]);
    end
else
    mov=zeros(height,width,1,tduration);
    for i=1:tduration
        % grab movie frame
        F = mvi.getFrame(i-1).getBuffer;

        % reshape into proper image
        F = reshape(F,3,width,height);
%         sum(sum(F(1,:,:) == F(2,:,:))) == height*width % sanity check: see if movie is really grayscale
%         sum(sum(F(1,:,:) == F(3,:,:))) == height*width
        F = F(1,:,:);
        mov(:,:,1,i) = permute(F,[3 2 1]);
    end
end


%test:
% imagesc(mov(:,:,1,1)) % same orientation as on stimulus monitor

