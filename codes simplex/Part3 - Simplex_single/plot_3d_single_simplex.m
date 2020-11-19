function [] = plot_3d_single_simplex(kappa_mat,est_coeff,d,TRI,Vertices_of_simplex)
load('training_and_validation_data.mat')
N_pts=200;
alpha=linspace(-0.2,1.2,N_pts)';
beta=linspace(-0.47,0.47,N_pts)';

alpha;%x axis vector
beta;%y axis vector
z = nan(length(alpha), length(beta));
x = nan(length(alpha), length(beta));
y = nan(length(alpha), length(beta));

%jj is x axis index
%ii is y axis index
for jj = 1:length(alpha)
    for ii = 1:length(beta)
         [tid,bary]   = tsearchn(Vertices_of_simplex, TRI, [alpha(jj) beta(ii)]);
         check=isfinite(tid);
         if check==1
         i_index=1;
         Basis=(factorial(d) ./ (factorial(kappa_mat(:,1)') ...
                    .*factorial(kappa_mat(:,2)').*factorial(kappa_mat(:,3))')) ...
                    .*((bary(i_index,1).^(kappa_mat(:,1)')) .* (bary(i_index,2).^(kappa_mat(:,2)')) ...
                    .* (bary(i_index,3).^(kappa_mat(:,3)')));
         

         z(ii, jj) = Basis*est_coeff;
         x(ii, jj) = alpha(jj);
         y(ii, jj) = beta(ii);
         end
    end
end


hold all, grid on, box on;
xlabel('alpha [deg]'); ylabel('beta[deg]'); zlabel('Cm');
surf(x,y,z,'EdgeColor','none','HandleVisibility','off')
plot3(alfa_true_tr,beta_tr,Cm_tr,'.','DisplayName','training')
plot3(alfa_true_tst,beta_tst,Cm_tst,'.','DisplayName','testing')
legend
title(['single simplex on 3d plot of degree: ' num2str(d)])
end

