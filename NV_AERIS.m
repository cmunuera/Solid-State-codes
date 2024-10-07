%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       AERIS sensor simulations                          %
%                       Carlos Munuera Javaloy                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear

% -things to load: bz, Bz, Om, noscs, oscpoints, points, tau
load("AERIS_DD_150.mat")

% Constants
%--------------------------------------------------------------------------
ge = 28.024e9;    % Hz/T
gh = 42.577478e6; % Nuclear gyromagnetic ratio (Hz/Tesla)
hbar = 1.054571628e-34;  % Hbar
T = 300;
kB = 1.38064852e-23;   % Julios * Kelvin-1 
D = 2.87e9;

% Simulation
%--------------------------------------------------------------------------
dt = 1/(Om*oscpoints*noscs);
meas_points = oscpoints*noscs;
mid_points = round(meas_points/2);
quarter_points = round(mid_points/2);
reps = size(bz, 1);

sequence = [quarter_points, mid_points+quarter_points];

OmNV = 20e6;
tpi = 1/(2*OmNV);

padding = 0;

t_vec = (1:reps)*tau;

% Operators
%--------------------------------------------------------------------------
sx = [0 1; 1 0];
sy = [0 -1i; 1i 0];
sz = [1 0; 0 -1];
sphi = @(phi) sx*cos(phi)+sy*sin(phi);
rho = [1 1; 1 1]/2;

H0 = @(bz) bz*ge*sz/2;
Hx = OmNV*sx/2;

% Loop
%--------------------------------------------------------------------------
spectrum = zeros(1, reps);

for m = 1:reps
    U = eye(2);

    for k = 1:meas_points
        if ismember(k, sequence)
            U = expm(-1i*2*pi*(H0(bz(m, k))+Hx)*tpi)*expm(-1i*2*pi*H0(bz(m, k))*(dt-tpi))*U;
        else
            U = expm(-1i*2*pi*H0(bz(m, k))*dt)*U;
        end
    end

    spectrum(m) = real(trace(U*rho*U'*sy));
end


%% Fourier
%--------------------------------------------------------------------------
[fspectrum, frecs] = fourier([spectrum, zeros(1, padding)], tau);

%% Plots
%--------------------------------------------------------------------------
figure(1)

subplot(2, 1, 1);
hold on
plot(t_vec, spectrum, 'LineWidth', 2)

subplot(2, 1, 2);
hold on
plot(frecs*sqrt(1+2*cos(55*pi/180)^2)/3, real(fspectrum.*exp(-1i*(2*pi*frecs*tau+pi/2))), 'LineWidth', 2)
%xlim([0, 500])
%ylim([0, 7e-5])


% Auxiliar functions
%--------------------------------------------------------------------------
function [fspectrum, frecs] = fourier(X, dt)
    L = length(X);
    Y = fftshift(fft(X));
    
    frecs = ((-L/2):((L/2)-1))/(L*dt);
    fspectrum = Y/L;
end
