%% Data Processing 4 Dof Robot -- Merbouche Mouloud & Tso Célia 
clear
close all 
clc

load('Data-20251124/q.mat');   
load('Data-20251124/tau.mat'); 

q4   = q_meas(1:4, :);   
Tau4 = Tau(1:4, :);      

%% Conversion des unités
q_conv = q4;                            
q_conv([1 2 4], :) = deg2rad(q4([1 2 4], :));  
q_conv(3, :)       = q4(3, :) / 1000;          

fc = 10;          
fe = 1000;       
Wn = 2 * fc / fe;
[b, a] = butter(1, Wn); 

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


save("q_conv.mat", "q_conv");   
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
