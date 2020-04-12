%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this script used to be "BleachingRate_USC.m"
clc
clear 

% Set physical constants
SpeedOfLight = 3e8;
PlanksConstant = 6.626e-34;
WaveLengthRange = 778:-2:380; % in nm
OpsinPhotosensitivity = 9.6e-21; % m^2 
wavelengths = 400:1:800;
WavelengthsCorrected = wavelengths * 1e-9; % nm -> m

% UDT Information
cd /Volumes/All_Staff/lab/Experiments/Calibration/2013-07-19
load cal_udt_spectrum;

spot_area = 0.011^2 * pi;
cd /Volumes/All_Staff/lab/Experiments/Calibration/2013-07-19
load LED_spectrum
LED_spectrum = LED_spectrum ./ norm(LED_spectrum);
RGBMatrix = LED_spectrum;

% Rod spectral sensitivty fit with Baylor nomogram
WL = 800:-1:400;
WL = WL.*0.001;
a0 = -5.2734;
a1 = -87.403;
a2 = 1228.4;
a3 = -3346.3;
a4 = -5070.3;
a5 = 30881;
a6 = -31607;
lambdaMax = 491;
LogPhotonSensitivity = a0*(log10((1./WL).*(lambdaMax/561))).^0 + a1*(log10((1./WL)*(lambdaMax/561))).^1 + a2*(log10((1./WL)*(lambdaMax/561))).^2 +  a3*(log10((1./WL)*(lambdaMax/561))).^3 + a4*(log10((1./WL)*(lambdaMax/561))).^4 + a5*(log10((1./WL)*(lambdaMax/561))).^5 + a6*(log10((1./WL)*(lambdaMax/561))).^6;
RodPhotonSensitivity = 10.^LogPhotonSensitivity;

%% LED

% greg confirmed that for 20200229 exp we used 80/20 mirror, just like xiaoyang calibration 20190112
powerNDF0_5 = [0.51e-6, 54.5e-9, 3.37e-9, 0.405e-9, 0.046e-9, 0.006e-9];

% Powers = powerNDF0_5(end); % Watts: Power made with 
PowerMeasured = powerNDF0_5;
PowerScaled = [powerNDF0_5(1), powerNDF0_5(1)/1e1, powerNDF0_5(1)/1e2, powerNDF0_5(1)/1e3, powerNDF0_5(1)/1e4, powerNDF0_5(1)/1e5];

for i = 1 : length(powerNDF0_5)
    Powers = PowerMeasured(i);
%     Powers = PowerScaled(i);
    TruePowerScalers = Powers ./ (cal_udt_spectrum(:,2)' * RGBMatrix');
    CalibratedRGBMatrix = TruePowerScalers .* RGBMatrix;
    SummedMonitorSpectra = CalibratedRGBMatrix';

    Powers
    CheckPower = dot(CalibratedRGBMatrix, cal_udt_spectrum(:,2))

    Intensity = (SummedMonitorSpectra .* WavelengthsCorrected') ./ (PlanksConstant * SpeedOfLight);
    PhotonFlux = Intensity ./ spot_area;
    AbsorptionRate = OpsinPhotosensitivity * dot(PhotonFlux,RodPhotonSensitivity([401:-1:1])');
    RodCollectingArea = 0.5e-12; % MOUSE m^2
    EffectivePhotonFlux = dot(PhotonFlux,RodPhotonSensitivity([401:-1:1])');
    PhotonCatchRate = EffectivePhotonFlux * RodCollectingArea;
    log_activated_rhodopsin_per_rod_per_sec(i) = log(PhotonCatchRate);
end