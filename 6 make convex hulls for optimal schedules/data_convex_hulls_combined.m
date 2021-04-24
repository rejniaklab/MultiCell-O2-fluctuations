function data_convex_hulls
  
  norm_data=0.2;
 
  for iter=1:4    
    switch iter
      case 1 
        file_nameI='data_black_influx'; 
        file_nameU='data_black_uptake'; 
        roi='black';
      case 2 
        file_nameI='data_blue_influx'; 
        file_nameU='data_blue_uptake'; 
        roi='blue';
      case 3 
        file_nameI='data_red_influx'; 
        file_nameU='data_red_uptake'; 
        roi='red';
      case 4 
        file_nameI='data_magenta_influx'; 
        file_nameU='data_magenta_uptake';    
        roi='magenta';
    end
        
    colH_U=[1,1,1]*0.5*0;   colL_U=[1,1,1]*0.95;
    colH_I=[0,153,0]/255; colL_I=[0,102,0]/255; 
    

    dataI=load(['results/',file_nameI,'.txt']);
    if max(dataI(:,1))-min(dataI(:,1))<0.1
      dataI(:,1)=dataI(:,1)+0.25*(0.5-rand(size(dataI,1),1));
    end
    dataU=load(['results/',file_nameU,'.txt']);
    if max(dataU(:,1))-min(dataU(:,1))<0.1
      dataU(:,1)=dataU(:,1)+0.25*(0.5-rand(size(dataU,1),1));
    end
    
    
    
    figure('position',[300,200,750,650])
    axis([0,100, 0,100, 0,5])
    hold on

    bin_all=[dataI(:,2),dataI(:,3),dataI(:,1)];
    [k_bin,vol]=convhulln(bin_all);
    trisurf(k_bin,dataI(:,2),dataI(:,3),dataI(:,1),'facecolor',...
            [153,255,255]/255,'edgecolor', 'none')
    alpha(0.25)
    
      
    ind=find(dataI(:,7)<norm_data); 
    if length(ind)>1
      if max(dataI(ind,1))-min(dataI(ind,1))<0.1
        dataI(ind,1)=dataI(ind,1)+0.25*(0.5-rand(length(ind),1));
      end
        
      bin_ind=[dataI(ind,2),dataI(ind,3),dataI(ind,1)];
      [k_bin,vol]=convhulln(bin_ind);
      trisurf(k_bin,dataI(ind,2),dataI(ind,3),dataI(ind,1),'facecolor',...
              colH_I,'edgecolor',colL_I)
      alpha(0.25)    
    end
  
  
    ind=find(dataU(:,7)<norm_data); 
    if length(ind)>1
      if max(dataU(ind,1))-min(dataU(ind,1))<0.1
        dataU(ind,1)=dataU(ind,1)+0.25*(0.5-rand(length(ind),1));
      end
        
      bin_ind=[dataU(ind,2),dataU(ind,3),dataU(ind,1)];
      [k_bin,vol]=convhulln(bin_ind);
      trisurf(k_bin,dataU(ind,2),dataU(ind,3),dataU(ind,1),'facecolor',...
              colH_U,'edgecolor',colL_U)
    end
    
    grid
    view(15,25)
    
    title(['Data: ',roi,':   cyan-all; 2-norm<',num2str(norm_data),...
           '; green: influx; black: uptake;'],'fontsize',15)
    xlabel('tumor fraction [%]','rotation',-3,'fontsize',18)
    ylabel('stromal fraction [%]','rotation', 50,'fontsize',18)
    zlabel('vascular fraction [%]','fontsize',18)
  
    print('-djpeg',['fig_',roi,'_',num2str(norm_data),'.jpg'])

  end
  
  
end
 