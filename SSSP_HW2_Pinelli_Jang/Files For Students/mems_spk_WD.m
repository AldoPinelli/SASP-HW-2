clear; close all; clc
%Hamin Jang - 10723759
%Aldo Pinelli - 10802164

%% Simscape File for Ground Truth

load('ssc_output.mat')

%% Sampling Frequency
fs = 192e3;

%% Sampling Period
Ts = 1/fs;

%% Simulation Duration
stop_time = 2;  % [seconds]

%% Input Signal

% Time Axis
t = 0:Ts:stop_time;
t = t(1:end-1);

% Signal Amplitude
A = 10;
% vin = A * sweeptone(1.99, 0.01, fs, 'SweepFrequencyRange', [500 20000]);
load('vin.mat');

%% Circuit Parameters
% Piezoelectric MEMS loudspeaker in Free-field conditions

% Transformer Ratio
alpha = 3.7e-4;
Seff = 2.0e-5;

% Electrical Domain
Re = 4;
Cp = 2.4e-8;

% Mechanical Domain
Rm = 9.7e-3;
Mm = 1e-6;
Cm = 2.2e-3;

% Acoustic Domain
Cbc = 3.6e-13; 
Ltube1 = 1e2;
Ltube2 = 1e2;
Ctube = 6.5e-13;
Rac = 5e6;

%% Removing Ideal Transformers (Mechanical Domain only)

gamma1 = alpha^-1;
gamma2 = Seff;

% Resistive Elements
R1 = Re/gamma1^2;
R2 = Rm;
R3 = Rac*gamma2^2;


% Dynamic Elements
% Capacitors
C1 = Cp*gamma1^2;
C2 = Cm;
C3 = Cbc/gamma2^2;
C4 = Ctube/gamma2^2;


% Inductors
L1 = Mm;
L2 = Ltube1*gamma2^2;
L3 = Ltube2*gamma2^2;

% Source
Fin = alpha * vin;

%% Setting of Free Parameters (Adaptation Conditions)
Z2 = R1;
Z5 = Ts/(2*C1); 
Z8 = R2;
Z9 = (2*L1)/Ts;
Z10 = Ts/(2*C2);
Z13 = Ts/(2*C3);
Z16 = (2*L2)/Ts;
Z19 = Ts/(2*C4);
Z22 = (2*L3)/Ts;
Z23 = R3;

Z21 = Z22 + Z23;
Z20 = Z21;
Z18 = (Z19*Z20)/(Z19+Z20);
Z17 = Z18;
Z15 = (Z17+Z16);
Z14 = Z15;
Z12 = (Z14*Z13)/(Z14+Z13);
Z11 = Z12;
Z7 = Z11 + Z8 + Z9 + Z10;
Z6 = Z7;
Z4 = (Z6*Z5)/(Z5+Z6);
Z3 = Z4;
Z1 = Z3 + Z2;




%% Computing Scattering Matrices
Bser = [1,1,1];
Qpar = [1,1,1];
Bser5 = [1,1,1,1,1];

Zser1= diag([Z1,Z2,Z3]);
Zpar1= diag([Z4,Z5,Z6]);
Zser2= diag([Z7,Z8,Z9,Z10,Z11]);
Zpar2=diag([Z12,Z13,Z14]);
Zser3= diag([Z15,Z16,Z17]);
Zpar3=diag([Z18,Z19,Z20]);
Zser4=diag([Z21,Z22,Z23]);

Sser1= eye(3) - 2*Zser1*Bser'*inv(Bser*Zser1*Bser')*Bser;
Sser2= eye(5) - 2*Zser2*Bser5'*inv(Bser5*Zser2*Bser5')*Bser5;
Sser3= eye(3) - 2*Zser3*Bser'*inv(Bser*Zser3*Bser')*Bser;
Sser4= eye(3) - 2*Zser4*Bser'*inv(Bser*Zser4*Bser')*Bser;

Spar1= 2*Qpar'*inv(Qpar*inv(Zpar1)*Qpar')*Qpar*inv(Zpar1) - eye(3);
Spar2= 2*Qpar'*inv(Qpar*inv(Zpar2)*Qpar')*Qpar*inv(Zpar2) - eye(3);
Spar3= 2*Qpar'*inv(Qpar*inv(Zpar3)*Qpar')*Qpar*inv(Zpar3) - eye(3);

disp(Sser1);


%% Initialization of Wavesa
a2=0;a8=0;a23=0; %Resistances
b5=0;b9=0;b10=0;b13=0;b16=0;b19=0;b22=0; % Inductor & Condeser


%% Initialization of Output Signals

Fout = zeros(1, length(t));

%% Simulation Algorithm

for n = 1 : length(Fin)

    %Manage Dynamic ELements
    a5 = b5;
    a10 = b10;
    a13 = b13;
    a19 = b19;

    a9 = -b9;
    a16 = -b16;
    a22 = -b22;
    a23=0;

    % Forward Scan
    b21 = Sser4(1,:) * [0; a22; a23];
    a20 = b21;

    b18 = Spar3(1,:) * [0;a19; a20];
    a17 = b18;

    b15 = Sser3(1,:) * [0; a16; a17];
    a14 = b15;

    b12 = Spar2(1,:) * [0; a13; a14];
    a11 = b12;

    b7 = Sser2(1,:) * [0; a8; a9; a10; a11];
    a6 = b7;

    b4 = Spar1(1,:) * [0; a5; a6];
    a3= b4;

    b1 = Sser1(1,:) * [0; a2; a3];

    % Local Root Scattering

    a1 = 2*Fin(n) - b1;

    % Backward Scan
    b3 = Sser1(3,:)*[a1; a2; a3];
    b2 = Sser1(2,:) * [a1; a2; a3];
    a4=b3;
    

    b5 = Spar1(2,:) * [a4; a5; a6];
    b6 = Spar1(3,:) * [a4; a5; a6];
    a7 = b6;
    

    b8 = Sser2(2,:) * [a7; a8; a9; a10; a11];
    b9 = Sser2(3,:) * [a7; a8; a9; a10; a11];
    b10 = Sser2(4,:) * [a7; a8; a9; a10; a11];
    b11 = Sser2(5,:) * [a7; a8; a9; a10; a11];
    a12 = b11;
    


    b13 = Spar2(2,:) * [a12; a13; a14];
    b14 = Spar2(3,:) * [a12; a13; a14];
    a15 = b14;
    

    b16 = Sser3(2,:) * [a15; a16; a17];
    b17 = Sser3(3,:) * [a15; a16; a17];
    a18 = b17;
    

    b20 = Spar3(3,:) * [a18; a19; a20];
    b19 = Spar3(2,:) * [a18; a19; a20];
    a21 = b20;
    

    b22 = Sser4(2,:) * [a21; a22; a23];
    b23 = Sser4(3,:) * [a21; a22; a23];
    
    

    % Read Output
    Fout(n) = (a23 + b23)/2;

end

%% Output Plots

% Computing acoustic pressure
pout = Fout ./ Seff;

% Time Domain Plots
figure
set(gcf, 'Color', 'w');
plot(t, pout, 'b', 'Linewidth', 2);
hold on
plot(gt(1, :), gt(2, :), 'r--', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$p_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Pressure - Time Domain','Fontsize',18,'interpreter','latex');

% Frequency domain Plots
nfft = 2^20;
res = fs/nfft;
f = (0:nfft/2-1) * res;

ir = impzest(vin/A, pout');
ir_gt = impzest(vin/A, gt(2, 1:end-1)');

tf = fft(ir, nfft);
tf_gt = fft(ir_gt, nfft);

abs_tf = abs(tf(1:nfft/2));
abs_tf_gt = abs(tf_gt(1:nfft/2));

figure
set(gcf, 'Color', 'w');
semilogx(f, 20*log10(abs_tf/2e-5), 'b', 'Linewidth', 2);
hold on
semilogx(f, 20*log10(abs_tf_gt/2e-5), 'r--', 'Linewidth', 2);
grid on;
xlim([500, 20000])
xlabel('Frequency [Hz]','Fontsize',16,'interpreter','latex');
ylabel('$\mathrm{SPL}\,[\mathrm{dB}_\mathrm{SPL}]$','Fontsize',16,'interpreter','latex');
legend('WDF','SSC','Fontsize',16,'interpreter','latex');
title('Output Sound Pressure Level - Frequency Domain','Fontsize',16,'interpreter','latex');

%% Error Plots

figure
set(gcf, 'Color', 'w');
hold on;
plot(t, pout - gt(2, 1:end-1), 'k', 'Linewidth', 2);
grid on;
xlabel('Time [seconds]','Fontsize',16,'interpreter','latex');
ylabel('$\mathcal{E}_{\mathrm{out}}$ [Pa]','Fontsize',16,'interpreter','latex');
title('Error Signal','Fontsize',16,'interpreter','latex');

%% Compute Mean Squared Error (MSE)

mse = mean((pout - gt(2, 1:end-1)).^2);
disp('MSE = ')
disp(mse)
