clc
clear
close all
dataname = '..\F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k')
N=length(Cm);%number of data points
beta = Z_k(:,2);  % measured angle of sideslip
load('..\Part2 - Kalman filter assignment\data_alfa_true.mat')
alfa_true=alfa_true';
PD = 0.7;  % percentage 70%
P=[alfa_true beta Cm];
idx = randperm(N);
Ptrain = P(idx(1:round(N*PD)),:);
Ptest = P(idx(round(N*PD)+1:end),:);
Cm_tr=Ptrain(:,3);
alfa_true_tr=Ptrain(:,1);
beta_tr=Ptrain(:,2);
Cm_tst=Ptest(:,3);
alfa_true_tst=Ptest(:,1);
beta_tst=Ptest(:,2);
N_tr=length(alfa_true_tr);
N_tst=length(alfa_true_tst);


save('training_and_validation_data.mat','N_tr','N_tst','N','Cm_tr','Cm_tst','alfa_true_tr','alfa_true_tst','beta_tr','beta_tst','alfa_true','beta')