spikes = datarun.spikes{cell_index};
NT = 216000;
Robs = NIM.Spks2Robs(spikes, binsize, NT );
XVi = XVi_old;

params_stim = NIM.create_stim_params([nLags NY, NX], 'stim_dt',	dt, 'upsampling', up_samp_fac, 'tent_spacing', tent_basis_spacing);
load mov_20180926_007.mat
[NY, NX, NFRAMES] = size(mov);
RF = datarun.vision.sta_fits{cell_index,1};
RF_center_x = RF.mean(1);
RF_center_y = NY - RF.mean(2);
re_size = 13/2 - 1; % reduced from 15 bc PC has low memory
x_upper = ceil(RF_center_x + re_size);
x_lower = floor(RF_center_x - re_size);
y_upper = ceil(RF_center_y + re_size);
y_lower = floor(RF_center_y - re_size);
mov_resize = mov(y_lower:y_upper, x_lower:x_upper, :);
NX = size(mov_resize,1); NY = NX;
mov = reshape(mov_resize, [NY*NX, NFRAMES])';
mov(mov < 0.5) = -0.48; % taken from xml contrast value
mov(mov > 0.5) = 0.48;
Xstim = NIM.create_time_embedding(mov, params_stim); % beware of memory

save mov_20180926_007_cell63.mat mov

[LLs,~,~,fp] = fit0.eval_model(Robs, Xstim, XVi );
LLs(2) = fitS.eval_model(Robs, Xstim, XVi );
LLs(3) = fit1.eval_model(Robs, Xstim, XVi );
LLs(4) = fit1S.eval_model(Robs, Xstim, XVi );
LLs(5) = fit2.eval_model(Robs, Xstim, XVi ); % fit2 (w nonlin) is not better than fit1
LLfit = LLs - fp.nullLL

[~,pred_rate0,~,~] = fit0.eval_model(Robs, Xstim, XVi );
[~,pred_rateS,~,~] = fitS.eval_model(Robs, Xstim, XVi );
[~,pred_rate1,~,~] = fit1.eval_model(Robs, Xstim, XVi );
[~,pred_rate1S,~,~] = fit1S.eval_model(Robs, Xstim, XVi );
[~,pred_rate2,~,~] = fit2.eval_model(Robs, Xstim, XVi );
pred_rates = [pred_rate0, pred_rateS, pred_rate1, pred_rate1S, pred_rate2];
Robs_test = Robs(XVi); 

% R2_part = 1 - mean( (spk_per_bin(bgn:fin)' - pred_rates(bgn:fin,4)) .^2) / var(spk_per_bin(bgn:fin))
for i = 1:5
    R2(i) = 1 - mean( (Robs_test - pred_rates(:,i)) .^2) / var(Robs_test);
end
R2
% plot(R2)
sum(R2>0) / size(R2,1) / size(R2,2) 

save('before_refit_63_13.mat', 'LLfit', 'R2', 'pred_rates', '-append')