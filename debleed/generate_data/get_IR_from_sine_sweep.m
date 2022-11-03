function ir = get_IR_from_sine_sweep(input, response, fs, nsamp)

    Nfft = 2^nextpow2(length(input) + length(response));
    ir = real(ifft(fft(response,Nfft)./fft(input,Nfft)));    
    ir = ir(1:nsamp);
    
    %plot spectrogram of IR
    figure;
    ftgram(ir,fs,'rir');
    title('Sine sweep IR');
     
 
end