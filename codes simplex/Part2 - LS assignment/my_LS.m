
%% load the data
clc
clear
close all
load('training_and_validation_data.mat')

vec_order=1:8;
RMS_training_error  =zeros(length(vec_order),1);
RMS_validation_error=zeros(length(vec_order),1);
VARmean             =zeros(length(vec_order),1);
coeff_matrix        =zeros(length(vec_order),sum(1:(length(vec_order)+1)));%store coeffcients at each order
validation_err_matrix=zeros(length(vec_order),N_tst);
for index_order=1:length(vec_order)
%start training
alfa=alfa_true_tr;
beta=beta_tr;
N=N_tr;
A=zeros(N_tr,sum(2:index_order));
for i_set=1:index_order
    dummy=zeros(N,i_set+1);
    for index=1:i_set+1
        dummy(:,index)=(alfa.^(i_set+1-index)).*(beta.^(index-1));
    end
    temp=size(dummy);
    start=i_set*(i_set+1)/2;
    A(:,start:start+temp(2)-1)=dummy;
end
clear dummy i_set index temp start alfa beta N
A=[ones(N_tr,1) A];
fi=(A'*A);
b=A'*Cm_tr;
my_coeff=fi\b;
length(my_coeff);% length of the coefficient vector
coeff_matrix(index_order,1:sum(1:(index_order+1)))=my_coeff;
Cm_eval=A*my_coeff;
training_err=Cm_eval-Cm_tr;
RMS_training_error(index_order)=rms(training_err)/(max(Cm_tr) - min(Cm_tr));
VARmean(index_order) = (mean(diag(inv(fi))));
clear A fi b

% evaluating the identified polynomial at validation datapoints
alfa=alfa_true_tst;
beta=beta_tst;
N=N_tst;
A=zeros(N,sum(2:index_order));
for i_set=1:index_order
    dummy=zeros(N,i_set+1);
    for index=1:i_set+1
        dummy(:,index)=(alfa.^(i_set+1-index)).*(beta.^(index-1));
    end
    temp=size(dummy);
    start=i_set*(i_set+1)/2;
    A(:,start:start+temp(2)-1)=dummy;
end
clear dummy i_set index temp start alfa beta N
A=[ones(N_tst,1) A];
Cm_eval=A*my_coeff;
validation_err=Cm_eval-Cm_tst;
validation_err_matrix(index_order,:)=validation_err;
RMS_validation_error(index_order)=rms(validation_err)/(max(Cm_tst) - min(Cm_tst));
[cr(:,:,index_order) , lags(:,:,index_order)] = xcorr(validation_err,100);
clear A b
end
% 95% Confidence Bounds
conf = 1.96/sqrt(length(Cm_tst));
clear validation_err
%% plot for polynomial surfaces from degree 1 to 6
figure
for order_select=1:6
my_coeff=coeff_matrix(order_select,1:sum(1:(order_select+1)));
my_coeff=my_coeff';
N_surf=100;
vec_alfa=linspace(-0.2,0.9,N_surf);
vec_beta=linspace(-0.2,0.25,N_surf);
[alfa_surf,beta_surf] = meshgrid(vec_alfa,vec_beta);
Cm_surf  = nan(length(vec_beta),length(vec_alfa));
for index_X=1:length(vec_alfa)
    alfa=alfa_surf(:,index_X);
    beta=beta_surf(:,index_X);
    N=N_surf;
    A=zeros(N,sum(2:order_select));
    for i_set=1:order_select
        dummy=zeros(N,i_set+1);
        for index=1:i_set+1
            dummy(:,index)=(alfa.^(i_set+1-index)).*(beta.^(index-1));
        end
        temp=size(dummy);
        start=i_set*(i_set+1)/2;
        A(:,start:start+temp(2)-1)=dummy;
    end
    clear dummy i_set index temp start alfa beta N
    A=[ones(N_surf,1) A];
    Cm_surf(:,index_X)=A*my_coeff;
    clear A b
end
subplot(3,2,order_select)
plot3(alfa_true_tr,beta_tr,Cm_tr,'.','DisplayName','training')
hold on
plot3(alfa_true_tst,beta_tst,Cm_tst,'.','DisplayName','testing')
surf(alfa_surf,beta_surf,Cm_surf,'EdgeColor','none','HandleVisibility','off')
xlabel('\alpha[rad]');ylabel('\beta[rad]');zlabel('C_m');grid on;
title(['polynomial of order:' num2str(order_select)])
end
%% plot Residual RMS 
figure
plot(vec_order,RMS_training_error,'DisplayName','training')
hold on
plot(vec_order,RMS_validation_error,'DisplayName','validation')
xticks(vec_order);legend;grid on;title('RMS relative error');xlabel('model order');ylabel('Rel RMS error')
%% plot validation error of order 3 and 6
val_order1=3;
val_order2=6;
figure;
subplot(1,2,1)
plot(validation_err_matrix(val_order1,:));title(['validation Model residual of order:' num2str(val_order1)]);grid on;xlabel('samples');ylabel('error')
subplot(1,2,2)
plot(validation_err_matrix(val_order2,:));title(['validation Model residual of order:' num2str(val_order2)]);grid on;xlabel('samples');ylabel('error')
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
%% Mean Parameter Variance
figure
plot(vec_order,log10(VARmean))
xticks(vec_order);grid on;title('Mean Parameter Variance');xlabel('model order');ylabel('log10 $$ (\overline {Var(parameters )} ) $$','Interpreter', 'Latex')
%% plot to compare coefficients using my OLS and MATLAB fit command on order 3
figure
order_compare=3;
sf = fit([alfa_true_tr, beta_tr],Cm_tr,['poly' num2str(order_compare) num2str(order_compare)]);
matalb_coeff=coeffvalues(sf);
my_coeff=coeff_matrix(order_compare,1:sum(1:(order_compare+1)));
error_bw_MATLAB_MyOLS=(my_coeff-matalb_coeff);
plot(error_bw_MATLAB_MyOLS)
grid on;title(['Difference between coefficients MyOLS and MATLAB polynomial fit for order:' num2str(order_compare)]);xlabel('coefficient index');ylabel('difference between coefficients')

