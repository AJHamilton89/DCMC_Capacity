% Calculates the Discrete-input Continuous-output Memoryless Channel (DCMC) capacity of AWGN and uncorrelated Rayleigh fading channels for BPSK, QPSK, 8PSK and 16QAM.
% Rob Maunder 18/05/2011
% Modifitcations Alex Hamilton 10/5/2021

% Copyright 2021 Alexander Hamilton, 2011 Robert G. Maunder. This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

clear all;

% +---------------------------------------------+
% | Choose the SNR, modulation and channel here |
% +---------------------------------------------+

% Affects the accuracy and duration of the simulation
symbol_count = 10000;

% Channel SNR
snr = -20:1:30; % dB


% Modulation scheme
% -----------------

% 2PSK
modulation = [+1, -1];

% 4PSK
% modulation = [+1, +i, -1, -i];

% 8PSK
% modulation = [+1, sqrt(1/2)*(+1+i), +i, sqrt(1/2)*(-1+i), -1, sqrt(1/2)*(-1-i), -i, sqrt(1/2)*(+1-i)];

% 16QAM
% modulation = sqrt(1/10)*[-3+3*i, -1+3*i, +1+3*i, +3+3*i, -3+1*i, -1+1*i, +1+1*i, +3+1*i, -3-1*i, -1-1*i, +1-1*i, +3-1*i, -3-3*i, -1-3*i, +1-3*i, +3-3*i];

% 32QAM
%modulation =       1/[(sqrt(1^2+1^2)+2*sqrt(1^2+3^2)+sqrt(3^2+3^2)+2*sqrt(5^2+3^2)+2*sqrt(5^2+1^2))/8]*...
                   [             -3+5*i, -1+5*i, 1+5*i, 3+5*i, ...
                         -5+3*i, -3+3*i, -1+3*i, 1+3*i, 3+3*i, 5+3*i, ...
                         -5+1*i, -3+1*i, -1+1*i, 1+1*i, 3+1*i, 5+1*i,...
                         -5-1*i, -3-1*i, -1-1*i, 1-1*i, 3-1*i, 5-1*i,...
                         -5-3*i, -3-3*i, -1-3*i, 1-3*i, 3-3*i, 5-3*i,...
                                 -3-5*i, -1-5*i, 1-5*i, 3-5*i];

% If you add more modulation schemes here, make sure their average transmit power is normalised to unity


% Channel
% -------

% Uncorrelated Rayleigh fading channel
%channel = sqrt(1/2)*(randn(1,symbol_count)+i*randn(1,symbol_count));

% AWGN channel
channel = ones(1,symbol_count);


% +------------------------+
% | Simulation starts here |
% +------------------------+

channel_capacity = zeros(size(snr));

C = zeros(size(snr));

for index = 1:length(snr)

    % Generate some random symbols
    symbols = ceil(length(modulation)*rand(1,symbol_count));

    % Generate the transmitted signal
    tx = modulation(symbols);

    % Generate some noise
    N0 = 1/(10^(snr(index)/10));
    noise = sqrt(N0/2)*(randn(1,symbol_count)+i*randn(1,symbol_count));

    % Generate the received signal
    rx = tx.*channel+noise;

    % Calculate the symbol probabilities
    probabilities = max(exp(-(abs(ones(length(modulation),1)*rx - modulation.'*channel).^2)/N0),realmin);

    % Normalise the symbol probabilities
    probabilities = probabilities ./ (ones(length(modulation),1)*sum(probabilities));

    % Calculate the capacity
    channel_capacity(index) = log2(length(modulation))+mean(sum(probabilities.*log2(probabilities)));

    C(index)=log2(1+(10^(snr(index)/10)));
    
end

%figure
plot(snr,channel_capacity);
xlabel('SNR [dB]');
ylabel('Capacity [bit/s/Hz]')

% figure
% plot(snr-10*log10(channel_capacity),channel_capacity);
% xlabel('E_b/N_0 [dB]');
% ylabel('Capacity [bit/s/Hz]')

