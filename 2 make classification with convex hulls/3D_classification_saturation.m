data=load('data_summary.txt');


 ind_1=find((data(:,4)>=0)&(data(:,4)<=12)); 
 
 ind_2=find((data(:,4)>12)&(data(:,4)<=24)); 
 
 ind_3=find((data(:,4)>24)&(data(:,4)<=36));  
 
 ind_4=find((data(:,4)>36)&(data(:,4)<=48));  

 ind_5=find((data(:,4)>48)&(data(:,4)<=60)); 
 
figure
axis([0,100, 0,100, 0,5])
hold on

  %plot3(data(ind_1,2),data(ind_1,3),data(ind_1,1),'bo', 'markerfacecolor','b')
  bin_data=[data(ind_1,2),data(ind_1,3),data(ind_1,1)];
  [k_bin,vol]=convhulln(bin_data);
  trisurf(k_bin,data(ind_1,2),data(ind_1,3),data(ind_1,1),'facecolor',...
         'b','edgecolor', 'none')
  alpha(1)


  %plot3(data(ind_2,2),data(ind_2,3),data(ind_2,1),'ro', 'markerfacecolor','r')
  bin_data=[data(ind_2,2),data(ind_2,3),data(ind_2,1)];
  [k_bin,vol]=convhulln(bin_data);
  trisurf(k_bin,data(ind_2,2),data(ind_2,3),data(ind_2,1),'facecolor',...
           'r', 'edgecolor', 'none')
  alpha(1)


  %plot3(data(ind_3,2),data(ind_3,3),data(ind_3,1),'yo', 'markerfacecolor','y')
  bin_data=[data(ind_3,2),data(ind_3,3),data(ind_3,1)];
  [k_bin,vol]=convhulln(bin_data);
  trisurf(k_bin,data(ind_3,2),data(ind_3,3),data(ind_3,1),'facecolor',...
         [1,128/255,0], 'edgecolor', 'none')
  alpha(1)

  
  %plot3(data(ind_3,2),data(ind_3,3),data(ind_3,1),'yo', 'markerfacecolor','y')
  bin_data=[data(ind_4,2),data(ind_4,3),data(ind_4,1)];
  [k_bin,vol]=convhulln(bin_data);
  trisurf(k_bin,data(ind_4,2),data(ind_4,3),data(ind_4,1),'facecolor',...
          'y', 'edgecolor', 'none')
  alpha(1)


  %plot3(data(ind_4,2),data(ind_4,3),data(ind_4,1),'o', 'markerfacecolor',[.5,.5,.5],'markeredgecolor',[.5,.5,.5] )
  bin_data=[data(ind_5,2),data(ind_5,3),data(ind_5,1)];
  [k_bin,vol]=convhulln(bin_data);
  trisurf(k_bin,data(ind_5,2),data(ind_5,3),data(ind_5,1),'facecolor',...
         [.85,.85,.85], 'edgecolor', 'none')
  alpha(1)

  
  axis([0,100, 0,100, 0,5])
  grid
  view(15,25)
  title({['Classification of tissue saturation:'],...
         ['[0-blue-12-red-24-orange-36-yellow-48-white-60] mmHg']},...
         'fontsize',15)
  xlabel('tumor fraction','rotation',0,'fontsize',15)
  ylabel('stromal fraction','rotation', 60,'fontsize',15)
  zlabel('vascular fraction','fontsize',15)
  
  print -djpeg 'tissue_saturation_classification'



