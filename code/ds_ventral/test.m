if scale == 1
    if order_post(t) <= ntrial(fid_seq(1))
        trial_flash_post = binned(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 :...
            trial_len*order_post(t)+section_idx(fid_seq(1),1));
    else
        trial_flash_post = binned(trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1) - trial_len/2 + 1 : ...
            trial_len*(order_post(t)-ntrial(fid_seq(1)))+section_idx(fid_seq(2),1));
    end
elseif scale == 2
    trial_flash_post = binned(trial_len*(trial_len*order_post(t)+section_idx(fid_seq(1),1) - trial_len/2 + 1 : ...
            trial_len*(trial_len*order_post(t)+section_idx(fid_seq(1),1));
end
other_flash_post = sum_flash_post - trial_flash_post;
mean_flash_post = other_flash_post ./ (length(trial_num_post) - 1) - mean_all .* length(trial_num_post) ./ (length(trial_num_post) - 1);

