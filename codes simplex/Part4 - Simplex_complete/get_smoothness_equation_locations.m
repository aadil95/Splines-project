function [locs] = get_smoothness_equation_locations(kappa_mat,kappa_m_mat)

locs=zeros(length(kappa_m_mat(:,1)),1);
index=1;
for i_kappa_mat=1:length(kappa_mat(:,1))
    for i_kappa_m_mat=1:length(kappa_m_mat(:,1))
        if kappa_mat(i_kappa_mat,:)==kappa_m_mat(i_kappa_m_mat,:)
        locs(index)=i_kappa_mat;
        index=index+1;
        end
    end
end    
    
clear index i_kappa_mat
end

