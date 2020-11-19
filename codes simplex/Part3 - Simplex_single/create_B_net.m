clc
clear
close all
n=2;
% d=4;
vec_d=[1 4 6 8];
figure
for index_d=1:4
d=vec_d(index_d);
kappa_mat=zeros(nchoosek(n+d,d),n+1);
index=1;
for i=d:-1:0
for j=d:-1:0
for k=d:-1:0
if(i+j+k)==d
kappa_mat(index,:)=[i j k];
index=index+1;
end
end
end
end
clear i j k index

c_plot=zeros(nchoosek(n+d,d),2);
index=1;
for i=1:d+1
    for j=1:i
        c_plot(index,:)=[i j];
    index=index+1;
    end
end
clear index i j
text_plot=char(kappa_mat+48);

subplot(2,2,index_d);
plot(c_plot(:,1),c_plot(:,2),'o', 'LineWidth',1,'MarkerFaceColor', 'b')
text(c_plot(:,1),c_plot(:,2),text_plot,'FontSize', 20);
% set(gca, 'color', 'none')
set(gca,'XTick',[])
set(gca,'YTick',[])
% set(gca,'visible','off')
title(['B-Net structure with degree:' num2str(d)])
end