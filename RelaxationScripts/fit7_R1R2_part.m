function []=fit7_R1R2_part(work_dir,fit_file,R1R2_file)

%% BASED on Eq. 1 and 2 in Farrow, Zhang, Szabo, Torchia and Kay, 1995
%% Note that no simplifying assumptions are needed since unlike in experiment,
%% all 5 frequencies are known (whereas the experiment measures 3 values only)
% Constants:
cd(work_dir);
u0 =  1.25663706E-6; % m kg s-2 A-2  [permeability of free space]
h = 6.62606957E-34; % m^2 kg /s  [h-planck]
Gamma_N = -4.3143721E+6; % Hz/T = Hz * A s^2 / kg = A s / kg [gyromagnetic ratio for N]
Gamma_H = 42.577130376E+6; % Hz/T = Hz * A s^2 / kg = A s / kg [gyromagnetic ratio for H]
r_H = 1.01E-10; % m  [distance between nuclei - assume constant]
Omega_H =  600.133E+6; % Hz;
Omega_N = 60.818E+6; % Hz;
PPM=160E-6;

% Load and Fourier Transform:
all=load(fit_file); 
N=size(all,1); % number of samples
N=N-mod(N,2); % make even
t_sec=all(1:N,1)*1E-9; % time in seconds
t_sec;
x=all(1:N,2); % fitted correlation
dT_sec = t_sec(2)-t_sec(1);
assert(abs(dT_sec-1E-12)<1e-15); % assume fitted functions in intervals of 1 ps

% Compute J spectral densities; see 5.95 in Protein NMR Spectroscopy: Principles and
% Practice, 2nd edition:
Fs=1/dT_sec; % sampling frequency, in Hz
xdft=real(fft(x))*dT_sec;  % in seconds; dT_sec accounts for d_tau  ;   NOTE: perhaps should be multiplied by 2pi? scaling wont affect R2/R1 ratio
J=xdft(1:N/2+1); % positive frequencies
J(2:end-1) = 2*J(2:end-1); % double all but zero and nyquist frequencies to account for the negative frequencies (that mirror the positive)
dF_Hz = Fs/N;
freq = 0:dF_Hz:Fs/2;

% Compute R1 and R2:
W0=1; % the first index in J is J(0)
WN=round(Omega_N/dF_Hz)+1; % normalized to freq space, +1 cause counting from 1
WH=round(Omega_H/dF_Hz)+1; % normalized to freq space, +1 cause counting from 1
d=(u0*h*Gamma_N*Gamma_H)/(8*pi*pi*r_H^3); % in Hz
c=(Omega_N/sqrt(3.0))*PPM; % in Hz
R1_DD= 0.25 *d^2*  (J(WH-WN) + 3*J(WN) + 6*(J(WH+WN))); % in Hz
R2_DD= 0.125*d^2*  (4*J(W0) + J(WH-WN) + 3*J(WN) + 6*J(WH) + 6*J(WH+WN)); % in Hz
R1_CSA=c^2*J(WN); % in Hz
R2_CSA=(c^2/6.0)*(3*J(WN)+4*J(W0)); % in Hz
R1=R1_DD+R1_CSA;
R2=R2_DD+R2_CSA; % ignore Ra as in Fawzi et al.

% Save:
R2onR1=R2/R1;
save(R1R2_file,'R1','R2','R2onR1','dF_Hz','d','c','-ascii');

end
