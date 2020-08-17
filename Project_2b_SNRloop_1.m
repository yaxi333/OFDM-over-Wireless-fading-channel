clear all; close all; clf

par.M = 4; % QPSK 
par.B = 1e6; % Bandwidth (tot)
par.fc = 2e9; % Carrier freq
par.lambda = 3e8/par.fc;
par.N0 = 2*2.07e-20; % Noise power [W/Hz](-164 dBm)
par.v = 15; % Relative velocity [m/s]
par.fD = par.v/par.lambda; % Doppler freq
par.Tcoh = 1/par.fD; % Coherence time
par.PL = 10^-10.1; % Path loss in linear (-101 dB)
par.tau = 5.4e-6; % Delay spread
par.P = 20:1:40; % Avg TX power [dBW] (per symbol)
par.Bcoh = 1/par.tau; % Coherence BW
par.O = 100; % Number of OFDM symbols
par.bps = log2(par.M); % Spectral eff / bits per symbol
par.channel = 'Fading_Channel'; % {'Fading_Channel' or 'awgn'}
par.Ts = 1/par.B; % Sample time
par.L = ceil(par.tau/par.Ts); % Minimum number of taps/Ncp

% Max number of samples, N + Ncp must be within Tcoh/10
par.N = floor(par.Tcoh/10/par.Ts); % Max no samples for a ~stationary channel
par.N = par.N - par.L; % Actual data samples w/o CP
par.N = pow2(floor(log2(par.N))); % actual data samples in power of 2 (# subcarriers)
par.R = par.N*par.bps/((par.N+par.L)*par.Ts); % Data rate


% Preallocate SER
SER = zeros(1,size(par.P,2));
SER_theory = zeros(1,size(par.P,2));
E_N0 = zeros(1,size(par.P,2));

for p = 1:length(par.P)
    % Create and scale symbols
    clear r
    clear h
    
    map = [ -1-1i, -1+1i, 1-1i, 1+1i ]; % Get the sombol map
    map = map./sqrt(mean(abs(map).^2)); % Normalize symbol energy, Es = 1
    bits = de2bi(0:par.M-1,par.bps,'left-msb'); % Create the bit pattern matching map
    E = 10^(par.P(p)/10)*par.Ts; % Energy per symbol
    map = map.*sqrt(E); % Scale to TX power

    % Create random data and map them to symbols
    data = randi(2,par.N*par.O,par.bps)-1;
    [~,idx] = ismember(data,bits,'rows'); % Find the corresponding symbol index
    S = map(idx); % SxO vector of symbols 

    % Serial to parallel
    S_p = reshape(S,[par.N par.O]); % Split up into subchannels 

    % Convert to time domain
    s_p = sqrt(par.N)*ifft(S_p,par.N,1); 

    % Add CP (precede s_p with end of itself)
    s = [s_p(par.N-par.L+1:par.N,:); s_p];

    % Parallel to serial, signal to transmit
    x = s.';

    tau = (0:par.L-1); %*round(par.tau/par.Ts);
    P = ones(1,par.L)*10^(par.P(p)/10)/par.L; % Power delay profile
    
    switch par.channel
        case 'Fading_Channel'

            % [r,h] = Fading_Channel(s, tau, fdts, P);
            % s = channel input (One OFDM symbol)
            % P = power delay profile [P/L]
            % tau = path delay in samples

            % Pass each OFDM symbol through the channel, r and h will be an OxS matrix
            for i = 1:par.O
                [r_temp,h_temp] = Fading_Channel(x(i,:), tau', par.fD*par.Ts, P);
                r(i,:) = r_temp;
                h(i,:) = h_temp(1,:);
            end

            % Add noise and path loss
            noise = sqrt(par.N0/2)*(randn(size(r))+1j*randn(size(r)));
            r = r*par.PL + noise;
            
            y = r.'; % Flip back
            h = h.'; % Flip back

            % discard trailing samples
            y = y(1:end-par.L+1,:);

            H = fft(h,par.N);

        case 'awgn'

            noise = sqrt(par.N0/2)*(randn(size(s))+1j*randn(size(s)));
            r = s*par.PL + noise;
            y = r;
            H = ones(par.N,par.O);
    end

    %--------------------------------------------------------------------------

    % Discard CP (The ISI guard samples)
    y = y(par.L+1:end,:);

    % Convert to freq domain
    Y = 1/sqrt(par.N)*fft(y);

    % Equalizer
    Y_hat = Y./H; 

    % Parellel to serial
    Y_s = reshape(Y_hat,1,[]);

    [~,idx] = min(abs(Y_s.' - map),[],2); % Get the closest points in the constellation
    S_hat = map(idx); % Estimated symbols

    SER(p) = numel(find(S_hat~=S))/numel(S_hat); % Symbol error rate 
    E_N0(p) = E*par.PL/par.N0/par.B;
    SER_theory(p) = 2*qfunc(sqrt(3*E_N0(p)/(par.M-1)));
    
end

scatterplot(Y_s); % Scatterplot
grid on
figure(1), semilogy(10*log10(E_N0),SER)
hold on, semilogy(10*log10(E_N0),SER_theory);
grid on
ylabel('SER')
xlabel('Es/N0 [dB]')
legend('Simulated QPSK','Theoretical QPSK');



