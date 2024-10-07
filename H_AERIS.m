%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Final sample simulation AERIS                        %
%                       Carlos Munuera Javaloy                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all

file_name = 'AERIS_DD_150.mat';

% Constants
%--------------------------------------------------------------------------
gh = 42.577478e6; % Nuclear gyromagnetic ratio (Hz/Tesla)
Bz = 2.1;
hbar = 1.054571628e-34;  % Hbar
mu0 = 4e-7*pi;  % N/A^2
T = 300;
kB = 1.38064852e-23;   % Julios * Kelvin-1 
therm = 2*pi*hbar*gh*Bz/(kB*T);

% Parameters
%--------------------------------------------------------------------------
% System
nuc_coords = ([-0.178/sqrt(3),      0,        0;
                0.178/(2*sqrt(3)),  0.178/2,  0;
                0.178/(2*sqrt(3)), -0.178/2,  0;
                0.178/sqrt(3),      0,        sqrt(0.248^2-0.178^2/3);
               -0.178/(2*sqrt(3)),  0.178/2,  sqrt(0.248^2-0.178^2/3)]*1e-9)';

N = size(nuc_coords, 2); % Number of nuclei
euler = [0, pi/5, pi/3]; % Rotate the molecule about ZYZ
CS = [3.66, 3.66, 3.66, 1.19, 1.19]; % Scalar chemical shifts

% Driving
Om = 150e3;
noscs = 2; % Oscillations per measurement stage
t_meas = noscs / Om; % Time it takes to perform a measuremnt

% Relaxation
T2 = 0.5;

% Simulation
t_tot = 0.5; % Total simulation time
oscpoints = 100; % Points per oscillation
dt = 1 / (oscpoints * Om);
points = 5000; 
meas_points = noscs * oscpoints;
mid_points = round(noscs*oscpoints/2);

tau = (t_tot - t_meas) / points; % Free times

% Operators
%--------------------------------------------------------------------------
[Ix, Iy, Iz, Id, Itotx, Itoty, Itotz] = operators(N);

rho0 = eye(2^N)/(2^N) + therm*Itotz/(2^N);

[Hcs, Hdd] = Hmol(nuc_coords, CS, euler, Ix, Iy, Iz, Bz, gh, mu0, hbar);
H0 = Hcs+Hdd;

Hx = Om * Itotx;
Hy = Om * Itoty;

% Master equation
ME = @(H) 2*pi*(-1i*kron(Id, H) + 1i*kron(H.', Id)) + 2/T2*kron(conj(Itotz), Itotz)...
         - 1/T2*kron(Id, Itotz'*Itotz) - 1/T2*kron(Itotz.'*conj(Itotz), Id);

% Superpropagators
%--------------------------------------------------------------------------
tpih = 1 / (4 * Om);

% Propagator of a Y pi/2 pulse
Uy = expm(ME(H0 + Hy)*tpih);

% Free propagator
U0 = expm(ME(H0)*tau);

% Drive propagator
Udpos = expm(ME(H0 + Hx)*dt);
Udneg = expm(ME(H0 - Hx)*dt);

%% AERIS Loop
%--------------------------------------------------------------------------
rhov = rho0(:);
rhov = Uy * rhov;

spectrumAE = zeros(points, meas_points);
for k = 1:points
    % Free time stage
    rhov = U0 * rhov;
    % Driving stage
    for ind = 1:mid_points
        rhov = Udpos * rhov;
        spectrumAE(k, ind) = real(trace(reshape(rhov, 2^N, 2^N)*Itotz));
    end
    for ind = mid_points+1:meas_points
        rhov = Udneg * rhov;
        spectrumAE(k, ind) = real(trace(reshape(rhov, 2^N, 2^N)*Itotz));
    end
end

hold on
plot(spectrumAE(5, :))
plot(spectrumAE(100, :))
plot(spectrumAE(139, :))

bz = spectrumAE * 1.3e-9*kB*T*2/(N*pi*hbar*gh);
save(file_name, "bz", "Bz", "Om", "noscs", "oscpoints", "points", "tau");