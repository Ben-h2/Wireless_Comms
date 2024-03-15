clear all
close all
clc

load("benchmark_parameter_174623_1472.mat")
load("benchmark_rece_data_174623_1472.mat")

k = 2048; %Total Number of Subcarriers
L = 200; %The number of zero padded symbols
Num_Taps = L+1; %Number of taps
Fc = 24E3;
B = 9E3; %Bandwidth
sampling_rate = 256E3;
W = 24; %The number of zero padded ofdm symbols in one packet
lambda = 24; %The over sampling factor

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passband filtering
Y_pb = bandpass(rece_data_ofdm_bench,[-1000+Fc,8000+Fc],sampling_rate);
figure()
plot(Y_pb)

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Doppler Rate a
v = 1.03; % Moving Speed
c = 1500; % Speed of Sound
a = v/c;
T_tx = 8.2695; % Transmitted signal duration
T_rx = 8.2687;
a_hat = T_tx/T_rx -1;

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resampling with a_hat

Y_pb_re = resample(Y_pb,round((1+a_hat)*1E5),1E5);

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Resampling from 256 kHz to 192 kHz

Ls = 192;
Ms = 256;
Lp = 24;
N = Lp*Ls -1 ;
h = Ls*fir1(N,1/Ms,kaiser(N+1,7.8562));
Y_pb_re_tilde = upfirdn(Y_pb_re,h,Ls,Ms);

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sychronization
load("pilot_signal_for_synchronization.mat")

StartPointCorr = xcorr(Y_pb_re_tilde,OFDM_data_pre_old);
figure()
plot(StartPointCorr);

[Max,n0] = max(StartPointCorr);% 2004510; %from xcorr plot

n0 = n0 - length(Y_pb_re_tilde);
figure()
subplot(2,1,1)
plot(Y_pb_re_tilde);
Y_pb_re_tilde = Y_pb_re_tilde(n0:end);
subplot(2,1,2)
plot(Y_pb_re_tilde);
title("Before and After Removing Padding")
%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 6 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Match Filtering P1

load("itc_1007_compfilter.mat")
Y_pb_re_mf_tilde =  conv(Y_pb_re_tilde,h_comp);