%% Fundamentals of Wireless Communication - Project
% Part 1 - task 1
% CP-OFDM
% Author: Tyler Bytendorp

% housekeeping
clear, clc
close all

%% Declare Constants

N = 16;                     % Number of subcarriers
OFDM_symbols = 4;           % Number of OFDM Symbols
Bits_per_QPSK_Symbol = 2;   % standard is 2 bits per QPSK symbol
h_td = [0.227; 0.46; 0.688; 0.46; 0.227];      % channel
L = 5;                      % Taps

%% Generate random bit data

b = rand(N.*Bits_per_QPSK_Symbol.*OFDM_symbols, 1);
a = find(b > 1/2);
c = find(b <= 1/2);
b(a) = 1;
b(c) = 0;
b_Tx = b;

%% convert bits to QPSK symbols

QPSK_syms_Tx = bits2QPSK(b_Tx);

%% Parse into column format
% Columns are each OFDM Symbol, Rows are each subcarrier QPSK symbol or 2 bits per subcarrier.

QPSK_syms_Tx2 = reshape(QPSK_syms_Tx, [N, OFDM_symbols]);
b_Tx2 = reshape(b_Tx, [N.*Bits_per_QPSK_Symbol, OFDM_symbols]);

%% Establish DFT/IDFT Matrix
n = 0:N - 1;        % overwritten later
[row, column] = meshgrid(n, n);
tt = row.*column;           % times table
F_N = exp(-1j.*2.*pi.*tt./N);
F_N = F_N./sqrt(N);

%% Start for loop to handle multiple OFDM symbols

QPSK_syms_Rx2 = zeros(size(QPSK_syms_Tx2));       % Pre-allocate space
for k = 1:OFDM_symbols

    %% Step 1: extract OFDM symbol
    
    d_fd_Tx = QPSK_syms_Tx2(:, k);
    
    %% Step 2: ifft(Frequency domain QPSK symbols) of Tx signal
    % Uses matrix method DFT - 1/sqrt(N) IDFT  - 1/sqrt(N);
    % can use fft and ifft
    
    %IFFT
    d_td = conj(F_N)*d_fd_Tx;
    
    %% Step 3: Add CP
    
    x = [d_td(end - (L - 1) + 1:end); d_td];
    
    %% Step 4: CP-OFDM Transmission
    
    h_sys = [flipud(h_td); zeros(N - (L - 1) - 1 ,1)];
    h_sys = circshift(h_sys, -1.*(L - 1), 1);
    
    y_td = zeros(N,1);
    for n = L:N + L - 1
         y_td(n - (L - 1)) = (h_sys.')*x(L:N + L - 1);
         h_sys = circshift(h_sys, 1, 1);
    end % end for
    
    %% Step 5: fft of Rx signal
    
    y_fd = F_N*y_td;
    
    %% Step 6: One-tap equalization
    
    % fft of h
    h_fd = vect_DFT(h_td, 16);
    
    d_fd_Rx = y_fd./h_fd;
    
    %% Reconstruct
    
    QPSK_syms_Rx2(:, k) = d_fd_Tx;

    %% end for loop to handle multiple OFDM symbols

end % end loop for multiple OFDM symbols

%% Reshape Rx QPSK symbols to vector format

QPSK_syms_Rx = reshape(QPSK_syms_Rx2, [N.*OFDM_symbols, 1]);

%% Convert QPSK symbols to bits

b_Rx = QPSK2bits(QPSK_syms_Rx);
b_Rx2 = reshape(b_Rx, [N.*Bits_per_QPSK_Symbol, OFDM_symbols]);

%% testing
aaa = [1+1j, 1-1j, -1+1j, -1-1j];
aaa = aaa.';
bbb = QPSK2bits(aaa);

%% Test/Check QPSK symbols

diff1 = QPSK_syms_Rx2 - QPSK_syms_Tx2;       % Should be 0's
disp(diff1)
disp('Difference in QPSK Symbols')

diff2 = b_Rx2 - b_Tx2;
disp(diff2)
disp('Difference in bits')


%% functions

%convert bits to QPSK symbols
function symbols=bits2QPSK(bits)

    if rem(length(bits),2)~=0
        bits=[bits, 0];
    end

    L=length(bits)/2;
    symbols=zeros(L,1);

    for k=1:L
        symbols(k)=-2*(bits(2*k)-0.5)+1j*-2*(bits(2*k-1)-0.5);
    end

    symbols = symbols/sqrt(2);

end

function X = vect_DFT(x, N)
    % x should be a column vector

    x = extendOrCutSignal(x, N);
    
    n = 0:N - 1;
    [kn, nk] = meshgrid(n, n);
    kn = kn.*nk;
    F_N = exp(-1j.*2.*pi.*kn./N);
    
    X = F_N*x;
    
    % below code sets angle to 0 if magnitude of dft is small.
    for k = 1:length(X)
        if abs(X(k)) < 1e-6
            X(k) = 0;
        end % end if
    end %end for

end % end function

function X = vect_IDFT(x, N)
    % x should be a column vector

    x = extendOrCutSignal(x, N);
    
    n = 0:N - 1;
    [kn, nk] = meshgrid(n, n);
    kn = kn.*nk;
    F_N = exp(-1j.*2.*pi.*kn./N);
    
    X = conj(F_N)*x;
    
    % below code sets angle to 0 if magnitude of dft is small.
    for k = 1:length(X)
        if abs(X(k)) < 1e-6
            X(k) = 0;
        end % end if
    end %end for

end % end function

function out = extendOrCutSignal(x, N)
    % must input column vector
    
    L = length(x);
    if N > L
        out = cat(1, x, zeros(N - L, 1));
    elseif N <= L
        out = x(1:end - (L - N));
    end % end if
end % end function

% convert QPSK symbols to bits
% NOTE: This function works for QPSK only and uses a slicer method.
function info_bits = QPSK2bits(info_symbols)
    
    % preallocate space
    L=length(info_symbols)*2;
    info_bits=zeros(L,1);        % column vector space allocated

    % Decode real and imaginary of symbols separately into bits
    I = real(info_symbols)<=0;
    Q = imag(info_symbols)<=0;
    
    % Parse bits back together
    for k=1:2:L - 1
        info_bits(k)=I((k+1)/2);
        info_bits(k+1)=Q((k+1)/2);
    end % end interleaving for loop
end % function