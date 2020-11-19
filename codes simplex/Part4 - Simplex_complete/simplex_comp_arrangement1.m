%% load the data
clc
clear
close all
load('training_and_validation_data.mat')
%% position of vertices annd creating the simplices 
V1=[min(alfa_true) min(beta)-0.14];
V2=[max(alfa_true)+0.4 0];
V3=[min(alfa_true) max(beta)+0.14];
V4=[0.25 0];
V=[V1;V2;V3;V4];
num_vertices=length(V);
vertices_names=arrayfun(@(n) {sprintf('V_{%d}', n)}, (1:num_vertices)');
T = delaunay(V);
T = sort(T,2);
TRI = triangulation(T,V);
n_simplices = size(TRI,1);
nb = neighbors(TRI);
%% plot training and testing and contain in simplices
figure
plot(alfa_true_tr,beta_tr,'.','DisplayName','training');hold all
plot(alfa_true_tst,beta_tst,'.','DisplayName','testing')
triplot(TRI,'r','HandleVisibility','off')
plot(V(:,1), V(:,2),'xk','LineWidth',2,'MarkerFaceColor', 'b','HandleVisibility','off')
text(V(:,1), V(:,2),vertices_names,'FontSize', 20);
legend;xlabel('alpha true[rad]');ylabel('beta[rad]')
title('simplex covering the entire data region (arrangement1)')
%plot the simplices at the centroids
IC = incenter(TRI);
simp_labels = arrayfun(@(v) {sprintf('t_{%d}', v)}, (1:n_simplices)');
text(IC(:,1),IC(:,2),simp_labels,'FontSize',20)
hold off
clear IC simp_labels
%% plot b coefficients on the simplices
d=4; %degree
for tj=1:n_simplices
[kappa_mat] = creat_multi_index_n_2(d);
bary_loc_coeff=kappa_mat./d;
ID=tj;
for i_coeff=1:length(bary_loc_coeff(:,1))
cart_loc_coeff{tj}(i_coeff,:) = barycentricToCartesian(TRI,ID,bary_loc_coeff(i_coeff,:));
coeff_names_for_loc{tj}(i_coeff,:)=['tj' num2str(tj) 'index ' char(kappa_mat(i_coeff,:)+48)];
end
end
for tj=1:n_simplices
figure; hold all
triplot(TRI,'r','HandleVisibility','off')
plot(V(:,1), V(:,2),'xk','LineWidth',2,'MarkerFaceColor', 'b','HandleVisibility','off')
text(V(:,1), V(:,2),vertices_names,'FontSize', 20);
xlabel('alpha true[rad]');ylabel('beta[rad]')
plot(cart_loc_coeff{tj}(:,1),cart_loc_coeff{tj}(:,2),'xr','LineWidth',3)
text(cart_loc_coeff{tj}(:,1),cart_loc_coeff{tj}(:,2),coeff_names_for_loc{tj},'FontSize',10,'HorizontalAlignment','right');
title(['plot B coefficient location in simplex tj:' num2str(tj) '(arrangement1)'])
end
clear tj d i_coeff ID coeff_names_for_loc cart_loc_coeff bary_loc_coeff kappa_mat
%% arrange the training and testing data with simplices
[t_tr,bary_tr]   = tsearchn(V, TRI, [alfa_true_tr beta_tr]);
[t_tst,bary_tst] = tsearchn(V, TRI, [alfa_true_tst beta_tst]);

%removing elements who are outside the simplices
t_tr=t_tr(all(isfinite(bary_tr(:,:)),2));
alfa_true_tr=alfa_true_tr(all(isfinite(bary_tr(:,:)),2));
beta_tr= beta_tr(all(isfinite(bary_tr(:,:)),2));
Cm_tr=Cm_tr(all(isfinite(bary_tr(:,:)),2));
bary_tr(all(isnan(bary_tr(:,:)),2),:) = [];
t_tst=t_tst(all(isfinite(bary_tst(:,:)),2));
alfa_true_tst=alfa_true_tst(all(isfinite(bary_tst(:,:)),2));
beta_tst= beta_tst(all(isfinite(bary_tst(:,:)),2));
Cm_tst=Cm_tst(all(isfinite(bary_tst(:,:)),2));
bary_tst(all(isnan(bary_tst(:,:)),2),:) = [];

%arranging them in order
Cmtr = [];
Cmtst = [];
alfa_truetst = [];
betatst = [];

for tj = 1:n_simplices
    b_tr{tj} = bary_tr(t_tr==tj,:);
    b_tst{tj}= bary_tst(t_tst==tj,:);
   
    % SORTING
    Cmtr = [Cmtr; Cm_tr(t_tr==tj)];
    Cmtst = [Cmtst; Cm_tst(t_tst==tj)];
    alfa_truetst = [alfa_truetst; alfa_true_tst(t_tst == tj)];
    betatst = [betatst; beta_tst(t_tst == tj)];
    Ntr{tj}=length(Cm_tr(t_tr==tj));
    Ntst{tj}=length(Cm_tst(t_tst==tj));
end
% clear tj Cm_tr Cm_tst alfa_true_tr alfa_true_tst beta_tr beta_tst d
clear tj d
%% estimating the coefficients
n=2;
vec_degree=[1 2 3 4 5 6 7 8 9 10 4 4 4];
vec_cont= [0 0 0 0 0 0 0 0 0  0 1 2 3];
for i_vec=1:length(vec_degree)
d=vec_degree(i_vec);%degree of each simplex
r=vec_cont(i_vec);%continuity order

[kappa_mat] = creat_multi_index_n_2(d);
total_bcoeffs_one_simp=nchoosek(d+n,n);
for tj=1:n_simplices
for i_index=1:Ntr{tj}
    B_tr{tj}(i_index,:)=(factorial(d) ./ (factorial(kappa_mat(:,1)') ...
                    .*factorial(kappa_mat(:,2)').*factorial(kappa_mat(:,3))')) ...
                    .*((b_tr{tj}(i_index,1).^(kappa_mat(:,1)')) .* (b_tr{tj}(i_index,2).^(kappa_mat(:,2)')) ...
                    .* (b_tr{tj}(i_index,3).^(kappa_mat(:,3)')));


end
clear i_index
for i_index=1:Ntst{tj}
    B_tst{tj}(i_index,:)=(factorial(d) ./ (factorial(kappa_mat(:,1)') ...
                .*factorial(kappa_mat(:,2)').*factorial(kappa_mat(:,3))')) ...
                .*((b_tst{tj}(i_index,1).^(kappa_mat(:,1)')) .* (b_tst{tj}(i_index,2).^(kappa_mat(:,2)')) ...
                .* (b_tst{tj}(i_index,3).^(kappa_mat(:,3)')));;

end
clear i_index
end
Btr = blkdiag(B_tr{:});
Btst = blkdiag(B_tst{:});
clear B_tr B_tst
%creating the smoothness matrix
H=[];
i_eqn=1;
% r=order of continuity
for m=0:r
   for tj=1:n_simplices
   for i_pair_row=1:length(nb(:,1)) %equations for each pair of simplices
   if (~isnan(nb(tj,i_pair_row)) && tj < nb(tj,i_pair_row))
   t1=nb(tj,i_pair_row);
   t2=tj;
   [non_edge_vertex_t1, loc_non_edge_t1]=setdiff(TRI(t1,:), TRI(t2,:));
   [non_edge_vertex_t2, loc_non_edge_t2]=setdiff(TRI(t2,:), TRI(t1,:));
   %create multi-indices for Left hand side
   [kappa_m_mat,k0,k1] = creat_multi_index_n_2_m(d,m,loc_non_edge_t1);
   %create multi-indices for Right hand side
   [kappa_zero_mat,gamma_mat] = creat_multi_index_n_2_zero(m,k0,k1,loc_non_edge_t2);  
   %coloumn locations of t1 in the kappa_mat matrix
   [t1_col_locs] = get_smoothness_equation_locations(kappa_mat,kappa_m_mat);

   bary_out_edge = cartesianToBarycentric(TRI, t2, V(non_edge_vertex_t1,:));%CHECK
   
   for i_row=1:length(t1_col_locs)
       H=[H; zeros(1,total_bcoeffs_one_simp*n_simplices)];
       H(i_eqn,total_bcoeffs_one_simp*(t1-1)+t1_col_locs(i_row))=-1;

       [t2_col_locs] = get_smoothness_equation_locations(kappa_mat,kappa_zero_mat{i_row});
        for rhs_row=1:length(t2_col_locs)
            H(i_eqn,total_bcoeffs_one_simp*(t2-1)+t2_col_locs(rhs_row))=+basis_func(m,gamma_mat(rhs_row,:),bary_out_edge);
        end
    i_eqn=i_eqn+1;
   end
   end
   end
   end
end
%arrange the smoothness matrix in the order of simplices
[rows_H, cols_H]=size(H);
F=Btr'*Btr;
M=[F H';H zeros(length(H(:,1)))];
if rank (H)< rows_H
    inv_M=pinv(M);
else
    inv_M=inv(M);
end
C1=inv_M(1:size(F),1:size(F));
var_param_log(i_vec)=log10(mean(diag(C1)));
est_coeff{i_vec}=C1*Btr'*Cmtr;% estimate the coeffients
%using the estimated coefficients on training data
Cm_est_tr=Btr*est_coeff{i_vec};
res_tr=Cm_est_tr-Cmtr;
residual_relative_rms_training(i_vec)=rms(Cm_est_tr-Cmtr)/(max(Cmtr)-min(Cmtr));
%using the estimated coefficients on validation data
Cm_est_tst=Btst*est_coeff{i_vec};
residual_validation{i_vec}=Cm_est_tst-Cmtst;
residual_relative_rms_validation(i_vec)=rms(Cm_est_tst-Cmtst)/(max(Cmtst)-min(Cmtst));
end
%% plot residuals as degree increases and continuity fixed
figure;hold all
x_plot=1:10';
y_plot_val=residual_relative_rms_validation(1:10)';
y_plot_tr =residual_relative_rms_training(1:10)';
plot(x_plot,y_plot_tr,'DisplayName','training');
plot(x_plot,y_plot_val,'DisplayName','validation');
xticks(x_plot);legend;grid on;title('RMS relative error with fixed continuity order:0 (arrangement1)');xlabel('model degree');ylabel('Rel RMS error')
%% Mean Parameter Variance as degree increases and continuity fixed
figure
x_plot=1:10';
y_plot=var_param_log(1:10);
plot(x_plot,y_plot)
xticks(x_plot);grid on;title('Mean Parameter Variance with fixed continuity order:0 (arrangement1)');xlabel('model degree');ylabel('log10 $$ (\overline {Var(parameters )} ) $$','Interpreter', 'Latex')
%% plot validation error of order 3 and 9
val_degree1=3;
val_degree2=9;
figure;
subplot(1,2,1)
plot(residual_validation{val_degree1});title(['validation Model residual of degree:' num2str(val_degree1) ' continuity order:0' '(arrangement1)']);grid on;xlabel('samples');ylabel('error')
subplot(1,2,2)
plot(residual_validation{val_degree2});title(['validation Model residual of degree:' num2str(val_degree2) ' continuity order:0' '(arrangement1)']);grid on;xlabel('samples');ylabel('error')
%% plot residuals as contunuity increases and degree fixed
select_index=[4 11 12 13];
vec_degree(select_index);
figure;hold all
x_plot=vec_cont(select_index);
y_plot_val=residual_relative_rms_validation(select_index)';
y_plot_tr =residual_relative_rms_training(select_index)';
plot(x_plot,y_plot_tr,'DisplayName','training');
plot(x_plot,y_plot_val,'DisplayName','validation');
xticks(x_plot);legend;grid on;title('RMS relative error with fixed degree: 4 (arrangement1)');xlabel('continuity order');ylabel('Rel RMS error')
%% plot the polynomial in 3d 
figure
for ii=1:length(select_index)
subplot(2,2,ii)
plot_3d_complete_simplex(est_coeff{select_index(ii)},vec_degree(select_index(ii)),vec_cont(select_index(ii)),TRI,V,n_simplices)
end
%% checking continuity order 0 of degree 4 simplex
coeff=est_coeff{4};
for tj=1:n_simplices
    simp_coeff_deg4{tj}=coeff((tj-1)*total_bcoeffs_one_simp+1:(tj)*total_bcoeffs_one_simp);
end
check_continuity=simp_coeff_deg4{2}-simp_coeff_deg4{3};