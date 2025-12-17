clear, close all, 
clc;
load validation.mat;

%% Filtrage 

q4va = q_meas_validation(1:4, :);   
Tau4val = tau_meas_validation(1:4, :);      

q_convval = q4val;                            
q_convval([1 2 4], :) = deg2rad(q4val([1 2 4], :));  
q_convval(3, :)       = q4val(3, :) / 1000;          

fc = 10;          
fe = 1000;       
Wn = 2 * fc / fe;
[b, a] = butter(1, Wn); 

[n_joints, n_samples] = size(q_convval);

qfval   = zeros(n_joints, n_samples);
dqfval  = zeros(n_joints, n_samples);
ddqfval = zeros(n_joints, n_samples);
taufval = zeros(n_joints, n_samples);

for i = 1:n_joints
    
    % Filtrage de la position 
    qfval(i, :) = filtfilt(b, a, q_convval(i, :));
    
    % Filtrage des couples
    taufval(i, :) = filtfilt(b, a, Tau4val(i, :));
    
    % Dérivation de la position filtrée
    dq_raw = derive(t, qfval(i, :));
    
    % Filtrage de la vitesse
    dqfval(i, :) = filtfilt(b, a, dq_raw);
    
    % Dérivation de la vitesse filtrée
    ddq_raw = derive(t, dqfval(i, :));
    
    % Filtrage de l'accélération
    ddqfval(i, :) = filtfilt(b, a, ddq_raw);
    
end

save("q_convval.mat", "q_convval");   
save("qfval.mat",     "qfval");
save("dqfval.mat",    "dqfval");
save("ddqfval.mat",   "ddqfval");
save("taufval.mat",   "taufval");


function df = derive(t, f)
    df = [f(1,2) - f(1,1), ...
        (f(1,3:end) - f(1,1:end-2)) / 2, ...
        f(1,end) - f(1,end-1) ...
        ] / (t(1,2) - t(1,1));
end
