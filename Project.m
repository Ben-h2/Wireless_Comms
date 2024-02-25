close all
clear all
clc

% Channel
h = [0.227 0.46 0.688 0.46 0.227];
L_taps = length(h);


% Declare number of subcarriers and number of OFDM symbols
K_subs = 16;
N_syms = 4;

% Create array of each symbol and its corresponding subcarriers
OFDM_syms = zeros(N_syms,K_subs);

% Fill the subcarriers with random constellation points
for n = 1:N_syms
    OFDM_syms(n,:) = createOFDMsymbols(K_subs);
    
    % figure()
    % plot(abs(OFDM_syms(n,:)))

    % disp(OFDM_syms(n,:))
    % figure()
    % scatter(real(OFDM_syms(n,:)), imag(OFDM_syms(n,:)))
end

% Convert to time domain
% Assuming a sample rate of 1Hz
OFDM_time = zeros(N_syms,K_subs);
for n = 1:N_syms
    OFDM_time(n,:) = ifft(OFDM_syms(n,:));
end

% Implement Cyclic Prefix
OFDM_cp = zeros(N_syms,K_subs+L_taps-1);
for n = 1:N_syms
    OFDM_cp(n,:) = [OFDM_time(n,K_subs-L_taps+2:K_subs) OFDM_time(n,:)];

    % figure()
    % plot(real(OFDM_cp(n,:)),'blue')
    % plot(imag(OFDM_cp(n,:)),'red')
end

% Combine all symbols into a single time domain wave form to be transmitted
OFDM_tx = OFDM_cp(1,:);
for n = 2:N_syms
    OFDM_tx = [OFDM_tx OFDM_cp(n,:)];
end

% Convolute with channel
OFDM_rx = conv(OFDM_tx,h);
for n = 1:N_syms
    OFDM_rx2(n,:) = conv(OFDM_cp(n,:),h);
end

% figure()
% plot(real(OFDM_rx),'blue')
% plot(imag(OFDM_rx),'red')

% Isolate the useful signal
OFDM_rx_useful = zeros(N_syms,K_subs);
OFDM_rx_useful2 = zeros(N_syms,K_subs);
for n = 1:N_syms
    p1 = n*(L_taps-1) + 1;
    p1_0 = (L_taps-1) + 1;
    p2 = p1+K_subs-1;
    p2_0 = p1_0+K_subs-1;
    OFDM_rx_useful(n,:) = OFDM_rx(p1:p2);
    OFDM_rx_useful2(n,:) = OFDM_rx2(n,p1_0:p2_0);

    % figure()
    % plot(real(OFDM_cp(n,:)),'blue')
    % plot(imag(OFDM_cp(n,:)),'red')
end

% Do fft
OFDM_rx_subs = zeros(N_syms,K_subs);
OFDM_rx_norm = zeros(N_syms,K_subs);

OFDM_rx_subs2 = zeros(N_syms,K_subs);
OFDM_rx_norm2 = zeros(N_syms,K_subs);

h_tilda = fft([h zeros(1,11)]);
for n = 1:N_syms
    OFDM_rx_subs(n,:) = fft(OFDM_rx_useful(n,:));
    OFDM_rx_norm(n,:) = OFDM_rx_subs(n,:) ./ h_tilda;

    OFDM_rx_subs2(n,:) = fft(OFDM_rx_useful2(n,:));
    OFDM_rx_norm2(n,:) = OFDM_rx_subs2(n,:) ./ h_tilda;

    % figure()
    % scatter(real(OFDM_rx_subs(n,:)), imag(OFDM_rx_subs(n,:)))

    % figure()
    % plot(abs(OFDM_rx_subs(n,:)))
    
    figure()
    plot(OFDM_rx_norm2(n,:))
end





