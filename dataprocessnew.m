clear
close all 
clc

load('Data-20251124/q.mat');   
load('Data-20251124/tau.mat'); 

q_rad = deg2rad(q_meas); 
fc = 10;        
fe = 1000;      
Wn = 2*fc/fe;   
[b, a] = butter(1, Wn); 
[n_joints, n_samples] = size(q_rad);

qf   = zeros(n_joints, n_samples);
dqf  = zeros(n_joints, n_samples);
ddqf = zeros(n_joints, n_samples);
tauf = zeros(n_joints, n_samples);
for i = 1:n_joints
    
    qf(i, :) = filtfilt(b, a, q_rad(i, :));
    
    tauf(i, :) = filtfilt(b, a, Tau(i, :));
    
    dq_raw = derive(t, qf(i, :));
    
    dqf(i, :) = filtfilt(b, a, dq_raw);
    
    ddq_raw = derive(t, dqf(i, :));
    
    ddqf(i, :) = filtfilt(b, a, ddq_raw);
    
end

save("q_rad.mat", "q_rad");
save("qf.mat", "qf");
save("dqf.mat", "dqf");
save("ddqf.mat", "ddqf");
save("tauf.mat", "tauf");

function df = derive(t,f)
df = [ f(1,2) - f(1,1) , (f(1,3:end) - f(1,1:end-2))/2 , f(1,end)-f(1,end-1)]/(t(1,2)-t(1,1));
end