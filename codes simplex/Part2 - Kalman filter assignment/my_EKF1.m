clc
clear
close all

%extract variables
dataname = '..\F16traindata_CMabV_2020';
load(dataname, 'Cm', 'Z_k', 'U_k')
dt = 0.01;%sample time
N=length(Cm);%number of data points
% measurements Z_k = Z(t) + v(t)
alpha_m = Z_k(:,1); % measured angle of attack
beta_m = Z_k(:,2);  % measured angle of sideslip
Vtot = Z_k(:,3);    % measured velocity

% input to Kalman filter
Au = U_k(:,1); % perfect accelerometer du/dt data
Av = U_k(:,2); % perfect accelerometer dv/dt data
Aw = U_k(:,3); % perfect accelerometer dw/dt data
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set initial values for states and statistics
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ex_0    = [150 0 0 0.3]'; % initial estimate of optimal value of x_k_1k_1
% Initial estimate for covariance matrix
P_0     = 5*diag([1 1 1 1]);

n       = 4; % number of states
nm      = 3; % number of measurements
m       = 3; % number of inputs
F       = [0 0 0 0;
           0 0 0 0;
           0 0 0 0;
           0 0 0 0]; %state matrix
B       = [1 0 0;0 1 0; 0 0 1;0 0 0]; % input matrix
G       = zeros(4,4); % noise input matrix

% System noise statistics:
Q = diag([1e-3 1e-3 1e-3 0]);

% Measurement noise statistics:
R = diag([0.035 0.013 0.110]);
%%
XX_k1k1 = zeros(n, N);
PP_k1k1 = zeros(n, N);
STDx_cor = zeros(n, N);
z_pred = zeros(nm, N);
IEKFitcount = zeros(N, 1);

x_k_1k_1 = Ex_0; % x(0|0)=E{x_0}
P_k_1k_1 = P_0; % P(0|0)=P(0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the Extended Kalman filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Extended Kalman Filter (EKF)
ti = 0; 
tf = dt;
% n = 1;%length(x_k); % n: state dimension
[dummy, Psi] = c2d(F, B, dt);   
[Phi, Gamma] = c2d(F, G, dt);   
% Run the filter through all N samples
for k = 1:N
    % Prediction x(k+1|k) 
    x_kk_1 = Phi*x_k_1k_1 + Psi*U_k(k,:)';

    % z(k+1|k) (predicted output)
    z_kk_1 = kf_calc_h(0, x_kk_1, U_k(k,:));
    z_pred(:,k) = z_kk_1;

    % Calc Phi(k+1,k) and Gamma(k+1, k)
%    precalculated as they are linear and same always
    
    % P(k+1|k) (prediction covariance matrix)
    P_kk_1 = Phi*P_k_1k_1*Phi' + Gamma*Q*Gamma'; 
    P_pred = diag(P_kk_1);
    stdx_pred = sqrt(diag(P_kk_1));

    % Correction
    Hx = kf_calc_Hx(0, x_kk_1, U_k(k,:)); % perturbation of h(x,u,t)
    % Pz(k+1|k) (covariance matrix of innovation)
    Ve = (Hx*P_kk_1 * Hx' + R); 
    % K(k+1) (gain)
    K = P_kk_1 * Hx' / Ve;
    % Calculate optimal state x(k+1|k+1) 
    x_k_1k_1 = x_kk_1 + K * (Z_k(k,:) - z_kk_1)'; 
    % P(k+1|k+1) (correction) using the numerically stable form of P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1; 
    P_k_1k_1 = (eye(n) - K*Hx) * P_kk_1 * (eye(n) - K*Hx)' + K*R*K';  
    P_cor = diag(P_k_1k_1);
    stdx_cor = sqrt(diag(P_k_1k_1));

    % Next step
    ti = tf; 
    tf = tf + dt;
    
    % store results
    XX_k1k1(:,k) = x_k_1k_1;
%     PP_k1k1(:,k) = diag(P_k_1k_1);
    STDx_cor(:,k) = stdx_cor;
end

figure
title('states');
subplot(1,4,1)
plot(XX_k1k1(1,:));
xlabel('samples')
title('reconstructed u')
subplot(1,4,2)
plot(XX_k1k1(2,:));
xlabel('samples')
title('reconstructed u')
subplot(1,4,3)
plot(XX_k1k1(3,:));
xlabel('samples')
title('reconstructed w')
subplot(1,4,4)
plot(XX_k1k1(4,:));
xlabel('samples')
title('reconstructed C_{alfa_{upwash}}')
C_alfaupwash_est=mean(XX_k1k1(4,5000:end));

figure
title('error of states')
subplot(1,4,1)
plot(STDx_cor(1,:));
xlabel('samples')
title('State estimation error u')
subplot(1,4,2)
plot(STDx_cor(2,:));
xlabel('samples')
title('State estimation error v')
subplot(1,4,3)
plot(STDx_cor(3,:));
xlabel('samples')
title('State estimation error w')
subplot(1,4,4)
plot(STDx_cor(4,:));
xlabel('samples')
title('State estimation error C_{alfa_{upwash}}')

figure
alfa_true=atan(XX_k1k1(3,:)./XX_k1k1(1,:));

title('alfa_m vs alfa_true')
plot(alpha_m,'Displayname','measured alfa');
hold on
plot(alfa_true,'Displayname','reconstructed true alfa');
ylabel('alfa [rad]')
xlabel('samples')
legend show
save('data_alfa_true.mat','alfa_true')