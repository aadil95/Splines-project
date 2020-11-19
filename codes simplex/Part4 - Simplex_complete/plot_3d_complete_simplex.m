function [] = plot_3d_complete_simplex(est_coeff,d,r,TRI,V,n_simplices)
[kappa_mat] = creat_multi_index_n_2(d);
total_bcoeffs_one_simp=nchoosek(d+2,2);
load('training_and_validation_data.mat')

for tj=1:n_simplices
    simp_coeff{tj}=est_coeff((tj-1)*total_bcoeffs_one_simp+1:(tj)*total_bcoeffs_one_simp);
end
N_pts=200;
alpha=linspace(-0.2,1.2,N_pts)';
beta=linspace(-0.4,0.4,N_pts)';
% [tid,bary]   = tsearchn(V, TRI, [alpha beta]);
% %removing elements who are outside the simplices
% tid=tid(all(isfinite(bary(:,:)),2));
% alpha=alpha(all(isfinite(bary(:,:)),2));
% beta= beta(all(isfinite(bary(:,:)),2));
% bary(all(isnan(bary(:,:)),2),:) = [];


alpha;%x axis vector
beta;%y axis vector
z = nan(length(alpha), length(beta));
x = nan(length(alpha), length(beta));
y = nan(length(alpha), length(beta));

%jj is x axis index
%ii is y axis index
for jj = 1:length(alpha)
    for ii = 1:length(beta)
         [tid,bary]   = tsearchn(V, TRI, [alpha(jj) beta(ii)]);
         check=isfinite(tid);
         if check==1
         i_index=1;
         Basis=(factorial(d) ./ (factorial(kappa_mat(:,1)') ...
                    .*factorial(kappa_mat(:,2)').*factorial(kappa_mat(:,3))')) ...
                    .*((bary(i_index,1).^(kappa_mat(:,1)')) .* (bary(i_index,2).^(kappa_mat(:,2)')) ...
                    .* (bary(i_index,3).^(kappa_mat(:,3)')));
         

         z(ii, jj) = Basis*simp_coeff{tid};
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

title(['simplex on 3d plot of degree: ' num2str(d) ' and continuity: ' num2str(r)])
end

