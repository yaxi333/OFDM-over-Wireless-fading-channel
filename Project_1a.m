 clear all, close all

% Parameter definition
par.fc = 2e9; % Carrier [GHz]
par.lambda = 3e8/par.fc; % Wavelength [m]
par.v = 30e3/3600; % Relative velocity [m/s]
par.Ts = 0.1e-3; % Sampling rate [s]
par.fD = par.v/par.lambda; % Doppler spread [Hz], no aliasing; BW=10k > 2fD=110
par.method = 'Spectrum'; % Choose method {'Filter', 'Spectrum'}
par.kc = 1; % Rician factor for dominating path, set to k = 0 for Rayleigh
par.Ns = 1e5; % Number of samples for the channel
par.N = par.Ns*2+1; % Number of samples for the filter g(t) 
par.fs = 1/par.Ts; % Sampling frequency

switch par.method
    
    case 'Filter'
        nTs = -par.N*par.Ts:par.Ts:par.N*par.Ts; % Time vector for filter

        % Define filter g(t)
        gt = besselj(1/4,2*pi*par.fD.*abs(nTs))./nthroot(abs(nTs),4); % For t~=0
        gt(par.N+1) = nthroot((pi*par.fD),4)/gamma(5/4); % For t=0
        gt = gt/sqrt(sum(abs(gt).^2)); % Normalise to unit energy
        Eg = sum(abs(gt).^2); % Check Energy of g(t)
        
        % Plot g(t) in time
        figure, plot(nTs.*1e3,gt);
        xlabel('Time [ms]')
        ylabel('g(t)')
        title('Time domain filter')
        
        % Generate Ns random complex distribution for the channel
        x = sqrt(1/2).*(randn(1,par.Ns)+1j*randn(1,par.Ns)); 
        var_x = var(x,[],2); % Check unit variance

        % Generate channel c(t)
        c = conv(x,gt,'same') + par.kc;
%         c = c./sqrt(mean(abs(c).^2)); % Renormalize to unit energy
        Ec = sum(abs(c).^2); % Check Energy of g(t)
        nTs = -par.Ns*par.Ts:par.Ts:par.Ns*par.Ts; % Time vector for the channel     

    case 'Spectrum'
        f = -par.fD:par.fs/par.Ns:par.fD; % Defining the frequency interval
        nTs = -par.Ns*par.Ts:par.Ts:par.Ns*par.Ts; % Time vector for sampling
        
        % Calculate Sf and Gf
        Sf = 1./(pi*par.fD.*sqrt(1-(f./par.fD).^2));
        Sf(1) = Sf(end); % Get rid of the 'Inf' at -fD
        Gf = sqrt(Sf);
        
        % The periodic extension of G(f) for k = 0:1:par.Ns-1
        Gfp = [Gf(length(Gf)/2:end) zeros(1,par.Ns-length(Gf)) Gf(1:length(Gf)/2-1)]; 
        var_Gfp = var(Gfp,[],2);
        
        % Random complex distribution
        a = sqrt(1/2/var_Gfp);
        X = a*(randn(1,par.Ns)+1j*randn(1,par.Ns)); % Random complex 
        var_X = var(X,[],2); % Check variance
        
        C = Gfp.*X;
        c = sqrt(par.Ns)*ifft(C,par.Ns,2);
     
        c = c + par.kc; % Add Rician k-factor
%         c = c./sqrt(mean(abs(c).^2,2)); % Renormalize to unit energy
        
end

% Plot channel distribution
figure, histogram(abs(c),20);
title('Distribution of c(t), Narrow band')

% Calculate expected values of c
E_Ec = mean(abs(c).^2,2); % E(energy) = 1
E_c = mean(real(c),2); % E(c) = 0 at Rayleigh

% Theoretical and estimated PSD
fD_vec = -par.fD:par.fs/par.Ns:par.fD;
Sc = 1./(pi*par.fD.*sqrt(1-(fD_vec./par.fD).^2));
Sc(1) = Sc(end);
% Sc_est = pwelch(c,[],[],[],par.fs); % Estimate of PSD
Sc_est = abs(sqrt(1/par.Ns)*fft(c)).^2/par.Ns;
figure
h1 = plot(fD_vec./par.fD,Sc,'LineWidth',2,'DisplayName','Theoretical');
grid on
hold on
h2 = plot((-par.Ns/2:par.Ns/2-1)*par.fs/par.Ns/par.fD,fftshift(Sc_est),'LineWidth',2,'DisplayName','Estimated');
axis([-1.1 1.1 0 max(get(h1,'Ydata'))])
title(['Doppler spectrum, ' par.method ' method'])
xlabel('f/f_D'), ylabel('S_c(f)')
grid on
legend

% PDF and CDF, theoretical and estimates of |c(t)|
x = 0:0.1:20; % x-axis def
dist = makedist('Rician','s',par.kc,'sigma',1); % Theoretical distribution
PDF = pdf(dist,x);
CDF = cdf(dist,x);
dist_est = fitdist(abs(c).','Rician'); % Estimated distribution
PDF_est = pdf(dist_est,x);
CDF_est = cdf(dist_est,x);
figure
subplot(1,2,1)
plot(x,PDF,'LineWidth',2,'DisplayName','Theory'), hold on
plot(x,PDF_est,'LineWidth',2,'DisplayName','Estimate')
title(['PDF k = ' num2str(par.kc)])
legend, grid on
subplot(1,2,2)
plot(x,CDF,'LineWidth',2,'DisplayName','Theory'), hold on
plot(x,CDF_est,'LineWidth',2,'DisplayName','Estimate')
title(['CDF k = ' num2str(par.kc)])
legend, grid on


% Theoretical and estimated ACF
Ac = besselj(0,2*pi*par.fD.*nTs) + par.kc^2; % Theoretical Ac(dt)
[Ac_est,lags] = xcorr(real(c),'unbiased'); % Estimated Ac(dt)
lags = lags(lags>=0);
Ac_est = Ac_est(lags>=0);
% Ac_est = Ac_est./max(Ac_est); % Normalize to max correlation = 1
figure
plot(nTs(par.Ns:end),Ac(par.Ns:end))
% title('ACF Theory')
% axis tight
% figure
hold on
plot(lags*par.Ts,Ac_est)
title(['ACF Estimate, k = ' num2str(par.kc)])
axis tight



