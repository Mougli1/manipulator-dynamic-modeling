clear
close all 
clc

load('Data-20251124/q.mat');   
load('Data-20251124/tau.mat');  % contient Tau


q4   = q_meas(1:4, :);   % positions des 4 premiers axes
Tau4 = Tau(1:4, :);      % couples des 4 premiers axes

%% Conversion des unités
% Joints 1, 2 et 4 : degrés -> radians
% Joint 3 : millimètres -> mètres

q_conv = q4;                            

q_conv([1 2 4], :) = deg2rad(q4([1 2 4], :));  %
q_conv(3, :)       = q4(3, :) / 1000;          % 3 en m

%% Paramètres du filtre passe-bas
fc = 10;          % fréquence de coupure (Hz)
fe = 1000;        % fréquence d'échantillonnage (Hz)
Wn = 2 * fc / fe; % fréquence normalisée
[b, a] = butter(1, Wn); 

%% Dimensions
[n_joints, n_samples] = size(q_conv);

%% Initialisation
qf   = zeros(n_joints, n_samples);
dqf  = zeros(n_joints, n_samples);
ddqf = zeros(n_joints, n_samples);
tauf = zeros(n_joints, n_samples);

for i = 1:n_joints
    
    % Filtrage de la position 
    qf(i, :) = filtfilt(b, a, q_conv(i, :));
    
    % Filtrage des couples
    tauf(i, :) = filtfilt(b, a, Tau4(i, :));
    
    % Dérivation de la position filtrée
    dq_raw = derive(t, qf(i, :));
    
    % Filtrage de la vitesse
    dqf(i, :) = filtfilt(b, a, dq_raw);
    
    % Dérivation de la vitesse filtrée
    ddq_raw = derive(t, dqf(i, :));
    
    % Filtrage de l'accélération
    ddqf(i, :) = filtfilt(b, a, ddq_raw);
    
end


save("q_conv.mat", "q_conv");   % q après conversion d'unités
save("qf.mat",     "qf");
save("dqf.mat",    "dqf");
save("ddqf.mat",   "ddqf");
save("tauf.mat",   "tauf");


function df = derive(t, f)
    df = [ ...
        f(1,2) - f(1,1), ...
        (f(1,3:end) - f(1,1:end-2)) / 2, ...
        f(1,end) - f(1,end-1) ...
        ] / (t(1,2) - t(1,1));
end
