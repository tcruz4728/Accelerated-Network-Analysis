% function wave = gen2PNwaveform_tau(fpos, ta, phase, fmin, fmax,tau0, tau1p5,
% 4datalen,initial_phase,snr,N, avec, normfac)
%Returns normalized phase vector of waveform in the Fourier Domain
function wave = gen2PNwaveform_tau(inParams,ta_index,phase,qcCoefs,snr)
fwavepos = waveform_tau(inParams,ta_index,phase,qcCoefs);

% fwavepos = waveform_tau(fpos,ta,phase,fmin,fmax,tau0,tau1p5,datalen,initial_phase, avec);

if mod(inParams.N,2) == 0
    fwaveneg = conj(fwavepos(end-1:-1:2));
else
    fwaveneg = conj(fwavepos(end:-1:2));
end

fwave = [fwavepos, fwaveneg];

wave = snr*inParams.normfac*fwave;




