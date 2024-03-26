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
T_rx = (1/sampling_rate)*2117317; % Ask TA how to calculate this value.
a_hat = 5E-5;%T_tx/T_rx -1;

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
Y_temp_save = Y_pb_re_tilde;
figure()
plot(StartPointCorr);

[~,n0xcorr] = max(StartPointCorr);% 2004510; %from xcorr plot

n0 = n0xcorr - length(Y_pb_re_tilde);
% n0 = 230000;
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
Y_pb_re_mf_tilde = Y_pb_re_mf_tilde(50+1:end);
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
Y_BB_bar = Y_BB_bar(2400+1:end); %Starting index found through inspection

%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid search

ninx = 0;
epsinx = 0;

for n_01 = 2200:2400
    ninx = ninx+1;
    for eps_1 = -2:0.1:2
        epsinx = epsinx+1;
        % 1-2 CFO compensation and Down Sampling
        n = n_01+[0:lambda:(k+L)*lambda-1];
        Y_BB_1_hat_down = Y_BB_bar(n).*exp(-1j*2*pi*eps_1*(n)*ts);

        % 3 Obtain the frequency domain 
        Z_1 = zp_fft(Y_BB_1_hat_down,k);

        % 4 calculate the power over null subcarriers
        P_null(ninx,epsinx)=0;
        for m = 1:length(k_null)
            P_null(ninx,epsinx)= P_null(ninx,epsinx) + abs(Z_1(k_null(m)))^2;
        end
    end
    epsinx = 0;
end
[~,inx] = min(P_null(:));
[ninx,epsinx]=ind2sub(size(P_null),inx);

n_start = 2200 + ninx;
disp(n_start)
eps_start = -2 + epsinx*0.1;
disp(eps_start)
figure();
h = pcolor(P_null');
set(h,'EdgeColor','None')
clear("P_null")
%%%%%%%%%%%%%%%%%%%%%%%%%%% STEP 10 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Grid Search for all subcarriers
n_start_W = n_start;
for W = 1:21
    for n_0w = n_start_W+(k+L)*lambda+[-2*lambda:2*lambda]
        ninx = ninx+1;
        for eps_1 = -2:0.1:2
            epsinx = epsinx+1;
            % 1-2 CFO compensation and Down Sampling
            n = n_0w+[0:lambda:(k+L)*lambda-1];
            Y_BB_1_hat_down = Y_BB_bar(n).*exp(-1j*2*pi*eps_1*(n)*ts);
    
            % 3 Obtain the frequency domain 
            Z_1 = zp_fft(Y_BB_1_hat_down,k);
    
            % 4 calculate the power over null subcarriers
            P_null(ninx,epsinx)=0;
            for m = 1:length(k_null)
                P_null(ninx,epsinx)=P_null(ninx,epsinx)+abs(Z_1(k_null(m)))^2;
            end
        end
        n_start_W = n_0w;
        epsinx = 0;
    end

[~,inx] = min(P_null(:));
[ninx,epsinx]=ind2sub(size(P_null),inx);

n_start_W = 2200 + ninx;
disp(ninx)
eps_start_W = -2 + epsinx*0.1;
disp(epsinx)
figure();
h = pcolor(P_null');
set(h,'EdgeColor','None')
clear("P_null")
end



% Timer
toc