function [kappa_zero_mat,gamma_mat] = creat_multi_index_n_2_zero(m,k0,k1,loc_non_edge_t2)

z=zeros(length(k0),1);
if loc_non_edge_t2==1
    kappa_mat_new(:,1)=z;
    kappa_mat_new(:,2)=k0;
    kappa_mat_new(:,3)=k1;
elseif loc_non_edge_t2==2
    kappa_mat_new(:,1)=k0;
    kappa_mat_new(:,2)=z;
    kappa_mat_new(:,3)=k1;
elseif loc_non_edge_t2==3
    kappa_mat_new(:,1)=k0;
    kappa_mat_new(:,2)=k1;
    kappa_mat_new(:,3)=z;
end        
gamma_mat = creat_multi_index_n_2(m);
for i=1:length(k0)
for i_g=1:length(gamma_mat(:,1))
    kappa_zero_mat{i}(i_g,:)=gamma_mat(i_g,:)+kappa_mat_new(i,:);
    end
end

clear i i_g
end

