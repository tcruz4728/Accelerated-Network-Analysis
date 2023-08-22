% function fwave = waveform_tau(fvec, t, phase, fmin, fmax,tau0, tau1p5,datalen,initial_phase,avec)
%Function to create Restricted 2PN Waveform in Fourier Domain
%Input is positive DFT frequency vector, returns phase of wave in fourier domain for
%postive DFT frequencies
function fwave = waveform(inParams, t,phase,gwCoefs)

fvec = inParams.fpos;
fmin = inParams.frange(1);
fmax = inParams.frange(2);
datalen = inParams.datalen;
initial_phase = inParams.initial_phase;
avec = inParams.avec;

% t = ta_index/inParams.Fs;

%% Constants
c = 3*10^8;
Msolar = 1.989*10^30;
G = 6.6743*10^-11;
cGfac = (c^3/G);
%% Calculate Mass Terms 
switch gwCoefs(3)
    case 1 %mass-space PSO
        m1 = gwCoefs(1)*Msolar;
        m2 = gwCoefs(2)*Msolar;
        M = (m1 + m2);
        u = m1*m2/M;
    case 2 %tau-space PSO
        tau0 = gwCoefs(1);
        tau1p5 = gwCoefs(2);
        M = (5/(32*fmin))*(tau1p5/(pi^2*tau0))*cGfac;
        u = (1/(16*fmin^2))*(5/(4*pi^4*tau0*tau1p5^2))^(1/3)*cGfac;
end
n = u/M;
%% Calculate Chirp Times
fminpifac = 2*pi*fmin;
gmfac = (M*pi*fmin*1/cGfac);
fminfac = (5/(96*fminpifac))*(1/n);

tau1 = fminfac*gmfac^(-1)*((743/336)+ (11*n/4));
% tau1 = (5/(192*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^-1)*(1/n)*((743/336)+ (11*n/4));
tau2 = 3/2*fminfac*gmfac^(-1/3)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));
% tau2 = (5/(128*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-1/3))*(1/n)*((3058673/1016064) + (5429*n/1008) + (617*n*n/144));
if gwCoefs(3) == 1
    tau0 = (5/(256*pi))*(1/fmin)*((G*M*pi*fmin/c^3)^(-5/3))*(1/n);
    tau1p5 = (1/8)*(1/fmin)*((G*M*pi*fmin/c^3)^(-2/3))*(1/n);
end

%% Alpha Terms
alpha0 = fminpifac*(3*tau0/5);

alpha1 = 0;

alpha2 = fminpifac*tau1;

alpha3 = -fminpifac*(3*tau1p5/2);

alpha4 = fminpifac*3*tau2;

alphaTerm = zeros(size(fvec));

alphaTerm(2:end) = alpha0*avec(1,:)... 
+ alpha1*avec(2,:)...
+ alpha2*avec(3,:)... 
+ alpha3*avec(4,:)... 
+ alpha4*avec(5,:);

F = 2*pi*fvec*(tau0 + tau1 - tau1p5 + tau2);

%% Final Phase Term
Psi = 2*pi*t*fvec - phase - pi/4 + alphaTerm + F + initial_phase;

%% Final Expression
fwave = exp(-1j*Psi);
% fwave = P;

%% Cut between fmin and fmax
min_index  = floor(datalen*fmin) + 1;
max_index = floor(datalen*fmax) + 1;

fwave(1:min_index-1) = 0;
fwave(max_index + 1: end) = 0;




