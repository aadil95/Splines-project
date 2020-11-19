%% plot for polynomial surfaces from degree 1 to 6
clc
clear
close all
load('my_estimated_single_simplex.mat')
n=2;
d=1;
my_coeff=B_coeff(1:nchoosek(n+d,d),d);
N_surf=1000;
vec_alfa=linspace(-0.1,1,N_surf)';
vec_beta=linspace(-0.2,1,N_surf)';
%convert to barry
[check_tr,barry_p]  = tsearchn(Vertices_of_simplex, simplex, [vec_alfa  vec_beta]);
% Selecting only that dataset which is inside the simplex
vec_alfa = vec_alfa(all(isfinite(barry_p(:,:)),2));
vec_beta = vec_beta(all(isfinite(barry_p(:,:)),2));
N_surf=length(vec_alfa);
barry_p(all(isnan(barry_p(:,:)),2),:) = [];
[alfa_surf,beta_surf] = meshgrid(vec_alfa,vec_beta);
Cm_surf  = nan(length(vec_beta),length(vec_alfa));

for xidx = 1 : length(vec_alfa)
  x = vec_alfa(xidx);
  for yidx = 1 : length(vec_beta)
    y = vec_beta(yidx);
    [check,barry]  = tsearchn(Vertices_of_simplex, simplex, [x  y]);
    z(xidx, yidx) = ((factorial(d) ./ (factorial(kappa_mat(:,1,d)') ...
                    .*factorial(kappa_mat(:,2,d)').*factorial(kappa_mat(:,3,d))')) ...
                    .*((barry(1).^(kappa_mat(:,1,d)')) .* (barry(2).^(kappa_mat(:,2,d)')) ...
                    .* (barry(3).^(kappa_mat(:,3,d)'))));
  return
  end
end
load('training_and_validation_data.mat')
figure
plot3(alfa_true_tr,beta_tr,Cm_tr,'o')
hold on
plot3(alfa_true_tst,beta_tst,Cm_tst,'o')
surf(alfa_surf,beta_surf,z)
xlabel('\alpha[rad]');ylabel('\beta[rad]');zlabel('C_m');grid on
title(['simplex of order:' num2str(d)])
