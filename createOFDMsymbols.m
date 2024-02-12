function [OFDM_symbols] = createOFDMsymbols(K)
% Create OFDM Signal from:
%   K --> Number of OFDM Subcarriers

% Generating random string of constellation points

% Generate two strings of random 0's and 1's
rand_phase1 = round(rand(1,K));
rand_phase2 = round(rand(1,K));
% Replace zeros with -1's
rand_phase1(rand_phase1==0)=-1;
rand_phase2(rand_phase2==0)=-1;

OFDM_symbols = 1/sqrt(2).*(rand_phase1 + 1j.*rand_phase2);
end

