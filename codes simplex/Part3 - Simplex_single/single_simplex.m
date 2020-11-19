
%% load the data
clc
clear
close all
load('training_and_validation_data.mat')
%% Part 2 position of vertices
y_offset_V1=-0.25;
y_offset_V2=+0.25;
x_offset_V3=+0.2;

% y_offset_V1=+0.1;
% y_offset_V2=-0.1;
% x_offset_V3=-0.5;
Vertices_of_simplex=[min(alfa_true) min(beta)+y_offset_V1; min(alfa_true) max(beta)+y_offset_V2;max(alfa_true)+x_offset_V3 0];

figure
plot(alfa_true_tr,beta_tr,'.','DisplayName','training');hold all
plot(alfa_true_tst,beta_tst,'.','DisplayName','testing')
plot(Vertices_of_simplex(:,1),Vertices_of_simplex(:,2),'xk','LineWidth',2,'MarkerFaceColor', 'b','DisplayName','Vertices of Simplex')
% simplex = delaunay(Vertices_of_simplex(:,1),Vertices_of_simplex(:,2));
simplex = delaunay(Vertices_of_simplex);
TRI = triangulation(simplex,Vertices_of_simplex);
triplot(simplex,Vertices_of_simplex(:,1),Vertices_of_simplex(:,2),'HandleVisibility','off');
legend;xlabel('alpha true[rad]');ylabel('beta[rad]')
title('simplex covering the entire data region')
%% Part 3 Identification of B coefficients

% convert the training and testing inputs to barrycentric coordinates
[check_tr,barry_tr]  = tsearchn(Vertices_of_simplex, simplex, [alfa_true_tr  beta_tr]);  
[check_tst,barry_tst] = tsearchn(Vertices_of_simplex, simplex, [alfa_true_tst beta_tst]);  

% Selecting only that dataset which is inside the simplex
alfa_true_tst = alfa_true_tst(all(isfinite(barry_tst(:,:)),2)); % alpha Simplex - Validation Data
beta_tst = beta_tst(all(isfinite(barry_tst(:,:)),2)); % beta Simplex - Validation Data
Cm_tr  = Cm_tr(all(isfinite(barry_tr(:,:)),2));     % Cm Simplex - Identification Data
Cm_tst = Cm_tst(all(isfinite(barry_tst(:,:)),2));	% Cm Simplex - Validation Data
N_tr=length(Cm_tr);
N_tst=length(Cm_tst);
barry_tr(all(isnan(barry_tr(:,:)),2),:) = [];
barry_tst(all(isnan(barry_tst(:,:)),2),:) = [];





% create matrix for indices
n=2;%simplex dimension
max_d=12;% degree
% B_tr=zeros(N_tr,nchoosek(n+max_d,max_d),max_d);
est_coeff=zeros(nchoosek(n+max_d,max_d),max_d);
log_mean_param_var=zeros(max_d,1);
res_tst=zeros(N_tst,max_d);
res_tr=zeros(N_tr,max_d);
for index_d=1:max_d

d=index_d;
% kappa_mat=zeros(nchoosek(n+d,d),n+1);
index=1;
for i=d:-1:0
for j=d:-1:0
for k=d:-1:0
if(i+j+k)==d
kappa_mat(index,:,d)=[i j k];
index=index+1;
end
end
end
end
clear i j k index

%training
% construct the regression matrix
B_tr=zeros(N_tr,nchoosek(n+d,d));
for i_index=1:N_tr
B_tr(i_index,:)=(factorial(d) ./ (factorial(kappa_mat(:,1,d)') ...
                    .*factorial(kappa_mat(:,2,d)').*factorial(kappa_mat(:,3,d))')) ...
                    .*((barry_tr(i_index,1).^(kappa_mat(:,1,d)')) .* (barry_tr(i_index,2).^(kappa_mat(:,2,d)')) ...
                    .* (barry_tr(i_index,3).^(kappa_mat(:,3,d)')));
end
% return

clear i_index
B_transp_B_inv=inv(B_tr'*B_tr);
log_mean_param_var(d)=log10(mean(diag(B_transp_B_inv)));
est_coeff(1:nchoosek(n+d,d),index_d)=B_transp_B_inv*B_tr'*Cm_tr;
Cm_pred_tr=B_tr*est_coeff(1:nchoosek(n+d,d),index_d);
res_tr(:,index_d)=Cm_pred_tr-Cm_tr;
RMS_training_error(index_d)=rms(Cm_pred_tr-Cm_tr)/(max(Cm_tr) - min(Cm_tr));

% clear B_tr

%validation
B_tst=zeros(N_tst,nchoosek(n+d,d));
for i_index=1:N_tst
%     B_tst(i_index,:)=factorial(d)*barry_tst(i_index,1).^(kappa_mat(:,1,d)').*barry_tr(i_index,2).^(kappa_mat(:,2,d)').*barry_tst(i_index,3).^(kappa_mat(:,3,d)')./((factorial(kappa_mat(:,1,d)').*factorial(kappa_mat(:,2,d)').*factorial(kappa_mat(:,3,d)')));    
      B_tst(i_index,:)=(factorial(d) ./ (factorial(kappa_mat(:,1,d)') ...
                    .*factorial(kappa_mat(:,2,d)').*factorial(kappa_mat(:,3,d))')) ...
                    .*((barry_tst(i_index,1).^(kappa_mat(:,1,d)')) .* (barry_tst(i_index,2).^(kappa_mat(:,2,d)')) ...
                    .* (barry_tst(i_index,3).^(kappa_mat(:,3,d)')));
end
clear i_index
Cm_pred_tst=B_tst*est_coeff(1:nchoosek(n+d,d),index_d);
res_tst(:,index_d)=Cm_pred_tst-Cm_tst;
RMS_validation_error(index_d)=rms(Cm_pred_tst-Cm_tst)/(max(Cm_tst) - min(Cm_tst));
[cr(:,:,index_d) , lags(:,:,index_d)] = xcorr(res_tst(:,index_d),100);

clear B_tst Cm_pred
end
conf = 1.96/sqrt(length(Cm_tst));
vec_degree=1:max_d;
%% plot Residual RMS 
figure
plot(vec_degree,RMS_training_error,'DisplayName','training')
hold on
plot(vec_degree,RMS_validation_error,'DisplayName','validation')
xticks(vec_degree);legend;grid on;title('RMS relative error');xlabel('model degree');ylabel('Rel RMS error')
%% plot validation error of order 3 and 9
val_degree1=3;
val_degree2=9;
figure;
subplot(1,2,1)
plot(res_tst(:,val_degree1));title(['validation Model residual of degree:' num2str(val_degree1)]);grid on;xlabel('samples');ylabel('error')
subplot(1,2,2)
plot(res_tst(:,val_degree2));title(['validation Model residual of degree:' num2str(val_degree2)]);grid on;xlabel('samples');ylabel('error')
%% Mean Parameter Variance
figure
plot(vec_degree,log_mean_param_var)
xticks(vec_degree);grid on;title('Mean Parameter Variance');xlabel('model degree');ylabel('log10 $$ (\overline {Var(parameters )} ) $$','Interpreter', 'Latex')
%% plot Auto correlation
val_order1=1;
val_order2=6;
figure;
subplot(1,2,1);hold all
plot(lags(:,:,val_order1),conf+lags(:,:,val_order1)*0,'--b','DisplayName','confidence interval')
plot(lags(:,:,val_order1),-conf+lags(:,:,val_order1)*0,'--b','HandleVisibility','off')
plot(lags(:,:,val_order1),cr(:,:,val_order1),'HandleVisibility','off');ylabel('Autocorrelation');xlabel('Number of Lags');title(['Autocorrelation of Model Residual of order:' num2str(val_order1)]);legend
subplot(1,2,2);hold all
plot(lags(:,:,val_order1),conf+lags(:,:,val_order2)*0,'--b','DisplayName','confidence interval')
plot(lags(:,:,val_order1),-conf-lags(:,:,val_order2)*0,'--b','HandleVisibility','off')
plot(lags(:,:,val_order2),cr(:,:,val_order2),'HandleVisibility','off');ylabel('Autocorrelation');xlabel('Number of Lags');title(['Autocorrelation of Model Residual of order:' num2str(val_order2)]);legend

%% plot 3d simplices
figure
for d=1:6
subplot(3,2,d)
plot_3d_single_simplex(kappa_mat(1:nchoosek(n+d,d),:,d),est_coeff(1:nchoosek(n+d,d),d),d,TRI,Vertices_of_simplex)
end
