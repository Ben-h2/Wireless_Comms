clear all
close all
clc
% Import Util
addpath('util');
addpath('data');

% Timer
tic

load("benchmark_parameter_174623_1472.mat")
load("benchmark_rece_data_174623_1472.mat")
load('ofdm_map.mat')

% Find Null SubK's
iszero = (ofdm_map == 0);
k_null = find(iszero);

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


%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Estimate Doppler Rate a
v = 1.03; % Moving Speed
c = 1500; % Speed of Sound
a = v/c;
T_tx = 8.2695; % Transmitted signal duration
T_rx = 8.2687; % Ask TA how to calculate this value.
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
figure();
plot(Y_pb_re_mf_tilde)

%removing 50 sample delay
Y_pb_re_mf_tilde = Y_pb_re_mf_tilde(50:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Passband to Baseband
fs = 192E3;
ts = 1/fs;

for n = 1:length(Y_pb_re_mf_tilde)
    Y_BB_I(n) = Y_pb_re_mf_tilde(n)*2*cos(2*pi*Fc*n*ts);
    Y_BB_Q(n) = Y_pb_re_mf_tilde(n)*2*sin(2*pi*Fc*n*ts);
end
Y_BB = Y_BB_I + Y_BB_Q.*1j;

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Match Filtering P2
beta = 0.125;
delay = 100;
span = 2*delay;
R = (rcosdesign(beta,span,lambda,"sqrt"));

Y_BB_bar = conv(Y_BB,R);

% plot(abs(Y_BB_bar))
Y_BB_bar = Y_BB_bar(4583:end); %Starting index found through inspection

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid search

ninx = 0;
epsinx = 0;

for n_01 = 2200:2400
    ninx = ninx+1;
    for eps_1 = -2:0.1:2
        epsinx = epsinx+1;
        % 1 CFO compensation
        for n = 1:(k+L)*lambda
            Y_BB_1_hat(n) = Y_BB_bar(n_01+n-1)*exp(-1j*2*pi*eps_1*(n_01+n-1)*ts);
        end
        % 2 Down-Sampling
        for i = 1:k+L
            Y_BB_1_hat_down(i) = Y_BB_1_hat(i*lambda);
        end
        % 3 Obtain the frequency domain 
        % for m = 1:k-1
        %     z(1,m) = 0;
        %     for i = 1:k+L
        %         Z_1(1,m) = z(1,m)+Y_BB_1_hat_down(i)*exp(-1j*2*pi*m*(i-1)/k);
        %     end
        % end
        Z_1 = zp_fft(Y_BB_1_hat_down,k);

        % 4 calculate the power over null subcarriers
        P_nuu(ninx,epsinx)=0;
        for m = 1:length(k_null)
            P_nuu(ninx,epsinx)= P_nuu(ninx,epsinx) + abs(Z_1(k_null(m)))^2;
        end
    end
    epsinx = 0;
end

% Timer
toc