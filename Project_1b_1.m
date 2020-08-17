% L no taps
clear all, close all 

% Parameter definition
par.fc = 2e9; % Carrier [GHz]
par.lambda = 3e8/par.fc; % Wavelength [m]
par.Ts = 0.1e-3; % Sampling rate [s]
par.Ns = 300; % Number of samples for the channel
par.fs = 1/par.Ts; % Sampling frequency
par.M = 64; % Zero-padding (oversampling in delay domain)
par.L = 3; % {1 2 3} Number of taps
par.kc = 0; % {0 1 10} Rician factor for dominating path, set k=0 for Rayleigh
par.fDTs = 0.005; % {0.1 0.005} For Clarke's spectrum, 0.4/fD = decorrelation time
par.fD = par.fDTs/par.Ts; % Doppler spread [Hz]

% Preallocation of c
cMN = zeros(par.M,par.Ns);

% ---------------Channel generation through Spectrum method----------------
f = -par.fD:par.fs/par.Ns:par.fD; % Defining the frequency interval
nTs = -par.Ns*par.Ts:par.Ts:par.Ns*par.Ts; % Time vector for sampling

% Generate G(f)
Sf = 1./(pi*par.fD.*sqrt(1-(f./par.fD).^2));
Sf = Sf(2:end-1); % Get rid of the 'Inf' at -fD and fD (not valid)
Gf = sqrt(Sf);

% The periodic extension of G(f) for k = 0:1:par.Ns-1
Gfp = [Gf(floor(length(Gf)/2):end) zeros(1,par.Ns-length(Gf)) Gf(1:floor(length(Gf)/2)-1)]; 
var_Gfp = var(Gfp,[],2);

% Random complex distribution
a = sqrt(1/2/var_Gfp); % Choose a such that the E[|c(t)|] = 1 ("eliminate the variance of G")
X = a*(randn(par.L,par.Ns)+1j*randn(par.L,par.Ns)); 
var_X = var(X,[],2); % Check unit variance

Cl = Gfp.*X; % Create channel in f domain
cl = sqrt(par.Ns)*ifft(Cl,[],2); % convert to t domain

cl = cl + par.kc; % Add Rician LOS factor
cl = cl./sqrt(mean(abs(cl).^2,2)); % Renormalize to unit expected energy
E_cl = mean(abs(cl).^2,2);

cMN(1:par.L,:) = cl; % Zero padding
CMN = fft(cMN,[],1); % FT in delay domain
% -------------------------------------------------------------------------
f = (0:par.M)./(par.Ts*par.Ns); % Frequency vector (oversampled)
tau = linspace(0,(par.L-1)*par.Ts/2,numel(Sf));
Ac_tau = sqrt(par.Ns)*ifft(Sf); % Autocorrelation in delay domain
DS_avg = sum(tau.*abs(Ac_tau).^2)/sum(abs(Ac_tau).^2); % delay spread in s

v = par.fD*par.lambda*3.6; % Relative velocity in kph

figure
subplot(1,2,1)
mesh(0:par.Ns-1,0:par.M-1,abs(CMN))
xlabel('Time [Ts]')
ylabel('Freq [1/(NsTs)]')
zlabel('|C(m,n)|')
title(['fDTs = ' num2str(par.fDTs) ', L = ' num2str(par.L) ', kc = ' num2str(par.kc)])
subplot(1,2,2)
surf(0:par.Ns-1,1:par.M,abs(CMN), 'MeshStyle', 'row'), view(2), axis tight
xlabel('Time [Ts]')
ylabel('Freq [1/(NsTs)]')
zlabel('|C(m,n)|')
title(['fDTs = ' num2str(par.fDTs) ', L = ' num2str(par.L) ', kc = ' num2str(par.kc)])

