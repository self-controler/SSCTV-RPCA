function HyperCube(data,band)
% % data_show=data(band,:,:);%%data的数据是band*length*width
% % [~,length,width]=size(data_show);
% % data_show_r=reshape(data_show(1,:,:),length,width);
% % data_show_g=reshape(data_show(2,:,:),length,width);
% % data_show_b=reshape(data_show(3,:,:),length,width);
% % data_show_r=(data_show_r-min(min(data_show_r)))/(max(max(data_show_r))-min(min(data_show_r)));
% % data_show_g=(data_show_g-min(min(data_show_g)))/(max(max(data_show_g))-min(min(data_show_g)));
% % data_show_b=(data_show_b-min(min(data_show_b)))/(max(max(data_show_b))-min(min(data_show_b)));
% % data_show=zeros(length,width,3);
% % data_show(:,:,1)=data_show_r;
% % data_show(:,:,2)=data_show_g;
% % data_show(:,:,3)=data_show_b;
% % figure;imshow(data_show);
data_show=double(data(:,:,band));%%data的数据是length*width*band
[length,width,~]=size(data_show);
data_show_r=reshape(data_show(:,:,1),length,width);
data_show_g=reshape(data_show(:,:,8),length,width);
data_show_b=reshape(data_show(:,:,17),length,width);
data_show_r=(data_show_r-min(min(data_show_r)))/(max(max(data_show_r))-min(min(data_show_r)));
data_show_g=(data_show_g-min(min(data_show_g)))/(max(max(data_show_g))-min(min(data_show_g)));
data_show_b=(data_show_b-min(min(data_show_b)))/(max(max(data_show_b))-min(min(data_show_b)));
data_show=zeros(length,width,3);
data_show(:,:,1)=data_show_r;
data_show(:,:,2)=data_show_g;
data_show(:,:,3)=data_show_b;
figure;imshow(data_show);
end