%% Initialization
% 13 clicks below highest setting reads 39.5 nW at 0.5 full screen, 0.1nW at 0
% full screen, 80nW at 1 full screen
% stage distance 7200um
% 50/50 reflecting mirror 4x objective
% rig A, lights off, curtain down
% Photodiode position was aligned to the objective axis at which the photometer
% reading was maximized by translating stage

my_path = '/Users/acquisition/Duke-devo';
%my_path = '/Users/acquisition/Photons';

addpath(genpath(my_path))
cd(my_path)

screen_number = 2; % Value = 0 (primary screen small), 1 (primary screen full), 2 (secondary screen full)
def_params = initialize_display('OLED', screen_number);


% set gamma OLED Aug 2, 2016
scale = [1.1399    1.0998    1.1027];
power = [1.1741    1.2998    1.3112];
offset = [-0.1445   -0.1023   -0.1054];
set_gamma_from_fit_params(scale, power, offset);
%% red

steps = 25;
run_gamma_flashes([scale, power, offset], steps, [1 0 0]);

r = [11.16 40.0 10.0 39.9 12.0 39.9 9.04 39.9 13.07 39.9 8.12 39.91 14.02 39.91 7.19 39.93 15.12 39.94 6.21 39.94 16.12 39.93 5.34 39.94 17.12 39.93 4.40 39.93 18.28 39.96 3.40 39.95 19.34 39.98 2.51 39.96 20.43 39.99 1.51 39.94 21.51 39.94 0.27 39.93 22.60 39.96 0.04 39.93 22.89 39.99];
%% green
run_gamma_flashes([scale, power, offset], steps, [0 1 0]);

g = [17.65 40.29 15.92 40.23 19.00 40.22 14.15 40.15 20.68 40.14 12.80 40.17 22.24 40.16 11.18 40.13 23.84 40.14 9.61 40.12 25.65 40.12 8.12 40.08 27.09 40.11 6.54 40.06 29.00 40.08 5.04 39.86 30.42 39.92 3.61 39.90 32.21 39.96 2.00 39.93 34.06 39.99 .33 39.95 35.60 39.99 0.04 39.95 36.62 39.99];
%% blue
run_gamma_flashes([scale, power, offset], steps, [0 0 1]);

b = [9.24 40.12 8.40 40.08 9.99 40.08 7.58 40.06 10.87 40.04 6.78 40.03 11.68 40.03 5.93 40.02 12.51 40.04 5.12 40.02 13.36 40.04 4.41 40.04 14.23 40.04 3.58 40.04 15.12 40.05  2.71 40.02 15.91 40.05 1.94 40.02 16.84 40.05 1.01 39.99 17.78 40.05 0.10 40.02 18.63 40.02 0.04 40.01 19.12 40.04];

%% fit
r = [11.16 40.0 10.0 39.9 12.0 39.9 9.04 39.9 13.07 39.9 8.12 39.91 14.02 39.91 7.19 39.93 15.12 39.94 6.21 39.94 16.12 39.93 5.34 39.94 17.12 39.93 4.40 39.93 18.28 39.96 3.40 39.95 19.34 39.98 2.51 39.96 20.43 39.99 1.51 39.94 21.51 39.94 0.27 39.93 22.60 39.96 0.04 39.93 22.89 39.99];
g = [17.65 40.29 15.92 40.23 19.00 40.22 14.15 40.15 20.68 40.14 12.80 40.17 22.24 40.16 11.18 40.13 23.84 40.14 9.61 40.12 25.65 40.12 8.12 40.08 27.09 40.11 6.54 40.06 29.00 40.08 5.04 39.86 30.42 39.92 3.61 39.90 32.21 39.96 2.00 39.93 34.06 39.99 .33 39.95 35.60 39.99 0.04 39.95 36.62 39.99];
b = [9.24 40.12 8.40 40.08 9.99 40.08 7.58 40.06 10.87 40.04 6.78 40.03 11.68 40.03 5.93 40.02 12.51 40.04 5.12 40.02 13.36 40.04 4.41 40.04 14.23 40.04 3.58 40.04 15.12 40.05  2.71 40.02 15.91 40.05 1.94 40.02 16.84 40.05 1.01 39.99 17.78 40.05 0.10 40.02 18.63 40.02 0.04 40.01 19.12 40.04];

%take out max values
real_r = r(1:2:end);
real_g = g(1:2:end);
real_b = b(1:2:end);

[~,i] = sort(get_sequence(steps));
sort_r = real_r(i);
sort_g = real_g(i);
sort_b = real_b(i);

[scale, power, offset] = fit_gamma([sort_r',sort_g',sort_b']);
scale = [1.03568066871491,1.03647726974954,1.03801579515696];
power = [1.07164946067619,1.05540770812658,1.03021347204221];
offset = [-0.00771159700138625,-0.0207938158112673,-0.0249832730848717];
%% compare to old values

%old values
scale_old = [1.1399    1.0998    1.1027];
power_old  = [1.1741    1.2998    1.3112];
offset_old  = [-0.1445   -0.1023   -0.1054];

figure;plot(scale_old,scale,'o');hold on
plot(power_old,power,'o');plot(offset_old,offset,'o');
refline(1)


xdata = 0:0.02:1;
for i = 1:3
    x1 = scale(i);
    x2 = power(i);
    x3 = offset(i);
    y_new = min(max(x1*xdata.^x2+x3,0),1);
    
    x1_old = scale_old(i);
    x2_old = power_old(i);
    x3_old = offset_old(i);
    y_old = min(max(x1_old*xdata.^x2_old+x3_old,0),1);
    
    figure;plot(xdata,y_new,'o-');hold on
    plot(xdata,y_old,'o-');
    title(num2str(i))
end



%%
