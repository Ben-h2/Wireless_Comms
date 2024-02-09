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