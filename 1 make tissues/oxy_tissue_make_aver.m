function [aver,seed]=oxy_tissue_make_aver(Vfrac,Tfrac,Sfrac)
% Vfrac -- fraction of tissue that is vasculature, i.e. 0.05 means 5%
% Tfrac -- fraction of tissue that is tumor, i.e. 0.2 means 20%
% Sfrac -- fraction of tissue that is stromal cells, i.e. 0.1 means 10%

aver=-1;

toSave=0;         % switch =1 will save data ==0 no data save
toSavelast=1;
Nsave=2500;       % frequency of saving the data
Ndraw=2500;       % frequency of drawing the data

pathDir=['data_V',num2str(Vfrac*100),'_T',num2str(Tfrac*100),'_S',...
         num2str(Sfrac*100)];

pathData=[pathDir,'/Data']; % name of the directiory for saving data
pathFigs=[pathDir,'/Figs']; % name of the directiory for saving figs    


% tissue not to dense
no_go=0;        
if (Vfrac+Tfrac+Sfrac>1) % if more than 100% do not simulate
  disp('more than 100% cellualrity & vascularity'); no_go=1;
end


% parameters
stiff=50;        % how hard the cells are pushing away from each other
nu=  250/10;     % viscosity of the medium
force_thr=5;     % threshold for stopping the forces
error_thr=1e-10; % threshold for stopping the oxygen stability

diamVess   =40;        % vessel diameter [microns]
diamCell   =15;        % tumor cell diameter [microns]
diamStroma = 7.5;      % stromal cell diameter

% domain 
xmin=-500; xmax=-xmin; ymin=xmin; ymax=-xmin;  % 1000x1000 micron^2=1mm^2
hg=5;           % grid width [microns]
dt=0.05;        % time step [sec]

% number of vessels and cells
Nvess=floor((xmax-xmin)*(ymax-ymin)*Vfrac/(pi*diamVess*diamVess/4));
Ntum=floor(((xmax-xmin)*(ymax-ymin))*Tfrac/(pi*diamCell*diamCell/4));
Nstroma=floor(((xmax-xmin)*(ymax-ymin))*Sfrac/(pi*diamStroma*diamStroma/4));
%Nvess,Ntum,Nstroma -- numbers of vessels, tumro cells and stromal cells

Doxy   =100;    % diffusion coefficient for oxygen [micron^2/s]
oxy_max= 60;    % level of oxygen in the vessel [mmHg=sigma-gram/micron^3]
oxy_init= 0.25; % initial oxygen level

oxy_upTum  =9.54/(hg*hg);   % max rate of oxy uptake by tumor cells
oxy_upStrom=9.54/(hg*hg);   % max rate of oxy uptake by stromal cells
km         =1.134;          % Michealis-Menten coefficient 


cfl_cond=Doxy*dt/(hg*hg);   % stability condition
if (cfl_cond>0.25)
  disp('stability condition is not satisfied'); pause;
end


% to remember the average oxygen and maximal force
avg_oxy=zeros(1,1);
max_for=zeros(1,1);
Nav=0;   % index for the above vectors

for ii=1:13 seed=floor(1000*rand); end
seed=floor(1000*rand); rng(seed); % rndom number generator



if no_go==0    % if the tissue is not too dense 
    
  % make output directory
  if (toSave==1)  mkdir(pathDir); mkdir(pathData); mkdir(pathFigs); end
  if (toSavelast==1)  mkdir(pathDir); mkdir(pathData); mkdir(pathFigs); end
  
  % random locations of vessels and cells   
  xyVess=zeros(Nvess, 2);
  vess_oxy=zeros(Nvess,1);
  for ii=1:Nvess
    xx=xmin+2*xmax*rand;
    yy=ymin+2*ymax*rand;
    xyVess(ii,1:2)=[xx,yy];
    vess_oxy(ii)=oxy_max;
  end

  xyTum=zeros(Ntum, 2);
  for ii=1:Ntum
    xx=xmin+2*xmax*rand;
    yy=ymin+2*ymax*rand;
    xyTum(ii,1:2)=[xx,yy];
  end

  xyStrom=zeros(Nstroma, 2);
  for ii=1:Nstroma
    xx=xmin+2*xmax*rand;
    yy=ymin+2*ymax*rand;
    xyStrom(ii,1:2)=[xx,yy];
  end


  % domain 
  [xx,yy]=meshgrid(xmin:hg:xmax,ymin:hg:ymax); % grid for oxygen 
  xxg=xx(1,:);  yyg=yy(:,1);  Ngx=length(xxg);  Ngy=length(yyg);


  % create the figure
  figure('position',[500,100,1500,500])
  pause(0.1)

  



% -------------- main code ---------------- %
  
  % initial oxygen level in the whole domain
  oxy=oxy_init*oxy_max*ones(Ngx+2,Ngy+2);

  
  iter=-1;                            % initial iteration  
 
  % save the previous level of oxygen to check convergence    
  old_oxy=10*oxy_max*ones(Ngx+2,Ngy+2);
  err_oxy=norm(old_oxy-oxy,2)/((Ngx+2)*(Ngy+2));  % change calculation 
  
  % conditions for cell relocation
  force_ok=0;                                       
  to_end=0; 
  
  
% ----- main loop ----- %  
 
  while (force_ok==0)||(err_oxy>error_thr) %force & error too large  
    iter=iter+1;   % increase iteration
           
    old_oxy=oxy;   % remember level of oxygen for convergence check
      
    % define influx of oxy from veins
    for ii=1:Nvess
      Nx=1+floor((xyVess(ii,1)-xmin)/hg);  % oxy grid closest to the vessel
      Ny=1+floor((xyVess(ii,2)-ymin)/hg);  
      for ix=-5:5                  % inspect near-by grid points
        for iy=-5:5                % to find this inside the vein
          if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)  
            ixy=sqrt((xyVess(ii,1)-(xmin+(Nx+ix)*hg))^2+...
                     (xyVess(ii,2)-(ymin+(Ny+iy)*hg))^2);
            if (ixy<0.5*diamVess)   % is inside the vein? 
              oxy(2+Nx+ix,2+Ny+iy)=vess_oxy(ii);
            end                    % influx rate
          end
        end
      end
    end
 
    
    % define diffusion of oxygen
    oxy(1:Ngx+2,1)      =oxy(1:Ngx+2,2);       % boundary conditions
    oxy(1:Ngx+2,Ngy+2)  =oxy(1:Ngx+2,Ngy+1);
    oxy(      1,1:Ngx+2)=oxy(      2,1:Ngx+2);
    oxy(  Ngx+2,1:Ngx+2)=oxy(  Ngx+1,1:Ngx+2);
    for ii=2:Ngx+1                             % diffusion loop 
      for jj=2:Ngy+1
        oxy(ii,jj)=oxy(ii,jj)+Doxy*(dt/hg^2)*(oxy(ii-1,jj)+oxy(ii+1,jj)+...
                   oxy(ii,jj-1)+oxy(ii,jj+1)-4*oxy(ii,jj));
      end
    end
    
    
    % define oxygen uptake by tumor cells
    for ii=1:Ntum
      Nx=1+floor((xyTum(ii,1)-xmin)/hg);  % grid closest to the tumor cell
      Ny=1+floor((xyTum(ii,2)-ymin)/hg);  
      for ix=-3:3                         % inspect near-by grid points
        for iy=-3:3                       % to find this inside the tumor
          if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)  
            ixy=sqrt((xyTum(ii,1)-(xmin+(Nx+ix)*hg))^2+...
                     (xyTum(ii,2)-(ymin+(Ny+iy)*hg))^2);
            if (ixy<0.5*diamCell)     
              %oxy(2+Nx+ix,2+Ny+iy)=(1-oxy_upTum*dt)*oxy(2+Nx+ix,2+Ny+iy);
              oxy(2+Nx+ix,2+Ny+iy)=max(0,oxy(2+Nx+ix,2+Ny+iy)*...
                              (1-dt*oxy_upTum/(km+oxy(2+Nx+ix,2+Ny+iy))));
            end                 % M-M method V*o/(km+o)
          end
        end   
      end
    end
    
    %oxy uptake by stromal cells
    for ii=1:Nstroma
      Nx=1+floor((xyStrom(ii,1)-xmin)/hg);  % grid closest to the cell
      Ny=1+floor((xyStrom(ii,2)-ymin)/hg);  
      for ix=-3:3                          % inspect near-by grid points
        for iy=-3:3                        % to find this inside the cell
          if (Nx+ix>0)&&(Nx+ix<=Ngx)&&(Ny+iy>0)&&(Ny+iy<=Ngy)  
            ixy=sqrt((xyStrom(ii,1)-(xmin+(Nx+ix)*hg))^2+...
                     (xyStrom(ii,2)-(ymin+(Ny+iy)*hg))^2);
            if (ixy<0.5*diamStroma)     
              oxy(2+Nx+ix,2+Ny+iy)=max(0,oxy(2+Nx+ix,2+Ny+iy)*...
                            (1-dt*oxy_upStrom/(km+oxy(2+Nx+ix,2+Ny+iy))));
            end             % uptake rate M-M method V*o/(km+o)
          end
        end   
      end   
    end
    
      
    % forces between Tum-Tum, Tum-Strom, Tum-Vess, Strom-Vess, Strom-Strom
    if (force_ok==0)  % previous force was too large to stop 
      forceTum=zeros(Ntum,2);    % all forces acting on tumor cells 
      forceStr=zeros(Nstroma,2); % all forces acting on stromal cells

      % Tumor <--> Tumor
      for ii=1:Ntum-1 %from first to last tum cell 
        for jj=ii+1:Ntum %inspects neighbor from current tum cell to last tum cell
          dx=xyTum(ii,1)-xyTum(jj,1);%positions between cells 
          dy=xyTum(ii,2)-xyTum(jj,2);%positons between cells
          dxy= sqrt(dx^2+dy^2);%distance between two cells, distance formula 
          if (dxy>0)&&(dxy<diamCell)% cant divide by 0 and dont want them to be too close
            forceTum(ii,1)=forceTum(ii,1)+stiff*(diamCell-dxy)*(dx/dxy);%changing force for ii nd jj tum cells
            forceTum(ii,2)=forceTum(ii,2)+stiff*(diamCell-dxy)*(dy/dxy);
            forceTum(jj,1)=forceTum(jj,1)-stiff*(diamCell-dxy)*(dx/dxy);% ii and jj are identical but opposite signs (same mag, diff direction)
            forceTum(jj,2)=forceTum(jj,2)-stiff*(diamCell-dxy)*(dy/dxy);
          end
        end   
      end
    
      % Tumor <--> Stroma
      restl=.5*diamCell+.5*diamStroma;
      for ii=1:Ntum %from first to last tum cell 
        for jj=1:Nstroma %inspects neighbor from current tum cell to last tum cell
          dx=xyTum(ii,1)-xyStrom(jj,1);%positions between cells 
          dy=xyTum(ii,2)-xyStrom(jj,2);%positons between cells
          dxy= sqrt(dx^2+dy^2);%distance between two cells, distance formula 
          if (dxy>0)&&(dxy<restl)% cant divide by 0 and dont want them to be too close
            forceTum(ii,1)=forceTum(ii,1)+stiff*(restl-dxy)*(dx/dxy);%changing force for ii nd jj tum cells
            forceTum(ii,2)=forceTum(ii,2)+stiff*(restl-dxy)*(dy/dxy);
            forceStr(jj,1)=forceStr(jj,1)-stiff*(restl-dxy)*(dx/dxy);%changing force for ii nd jj tum cells
            forceStr(jj,2)=forceStr(jj,2)-stiff*(restl-dxy)*(dy/dxy);        
          end
        end   
      end

      % Tumor <-- Vessel
      restl=.5*diamCell+.5*diamVess;
      for ii=1:Ntum %from first to last tum cell 
        for jj=1:Nvess %inspects neighbor from current tum cell to last tum cell
          dx=xyTum(ii,1)-xyVess(jj,1);%positions between cells 
          dy=xyTum(ii,2)-xyVess(jj,2);%positons between cells
          dxy= sqrt(dx^2+dy^2);%distance between two cells, distance formula 
          if (dxy>0)&&(dxy<restl)% cant divide by 0 and dont want them to be too close
            forceTum(ii,1)=forceTum(ii,1)+stiff*(restl-dxy)*(dx/dxy);%changing force for ii nd jj tum cells
            forceTum(ii,2)=forceTum(ii,2)+stiff*(restl-dxy)*(dy/dxy);
          end
        end   
      end

     % push tumor cells 
     xyTum=xyTum+forceTum*dt/nu;
    
    
    
     % Stroma <--> Stroma      
     for ii=1:Nstroma-1 %from first to last tum cell 
       for jj=ii+1:Nstroma %inspects neighbor from current tum cell to last tum cell
         dx=xyStrom(ii,1)-xyStrom(jj,1);%positions between cells 
         dy=xyStrom(ii,2)-xyStrom(jj,2);%positons between cells
         dxy= sqrt(dx^2+dy^2);%distance between two cells, distance formula 
         if (dxy>0)&&(dxy<diamStroma)% cant divide by 0 and dont want them to be too close
           forceStr(ii,1)=forceStr(ii,1)+stiff*(diamStroma-dxy)*(dx/dxy);%changing force for ii nd jj tum cells
           forceStr(ii,2)=forceStr(ii,2)+stiff*(diamStroma-dxy)*(dy/dxy);
           forceStr(jj,1)=forceStr(jj,1)-stiff*(diamStroma-dxy)*(dx/dxy);% ii and jj are identical but opposite signs (same mag, diff direction)
           forceStr(jj,2)=forceStr(jj,2)-stiff*(diamStroma-dxy)*(dy/dxy);
         end
       end   
     end
   
     % Stroma <-- Vessel      
     restl=.5*diamStroma+.5*diamVess;
     for ii=1:Nstroma %from first to last tum cell 
       for jj=1:Nvess %inspects neighbor from current tum cell to last tum cell
         dx=xyStrom(ii,1)-xyVess(jj,1);%positions between cells 
         dy=xyStrom(ii,2)-xyVess(jj,2);%positons between cells
         dxy= sqrt(dx^2+dy^2);%distance between two cells, distance formula 
         if (dxy>0)&&(dxy<restl)% cant divide by 0 and dont want them to be too close
           forceStr(ii,1)=forceStr(ii,1)+stiff*(restl-dxy)*(dx/dxy);%changing force for ii nd jj tum cells
           forceStr(ii,2)=forceStr(ii,2)+stiff*(restl-dxy)*(dy/dxy);
         end
       end   
     end
     
     % push atromal cells
     xyStrom=xyStrom+forceStr*dt/nu;
    
     
     % if forces were small enough force_ok=1 and forces will not be 
     % calculated in teh next iteration
     if force_ok==0
       maxF=max(max(max(forceStr)),max(max(forceTum)));  
       if (maxF<=force_thr) force_ok=1; end
     end 
     
      
     
% push back the vessels if they are oouside the domain
    for ii=1:Nvess    
      if (xyVess(ii,1)<xmin) xyVess(ii,1)=xmin+2*diamVess;end
      if (xyVess(ii,2)<ymin) xyVess(ii,2)=ymin+2*diamVess;end
      if (xyVess(ii,1)>xmax) xyVess(ii,1)=xmax-2*diamVess;end
      if (xyVess(ii,2)>ymax) xyVess(ii,2)=ymax-2*diamVess;end  
    end
   
    
  end  % if (force_ok==0) -- end of applying forces 
 
 
 
  % calculate the norm-2 difference between old and new oxy levels 
  % to check whether the difference is small (convergence) 
  err_oxy=norm(old_oxy-oxy,2)/(Ngx*Ngy);
    
 
  % if forces and error are small emough, to_end=1
  if (force_ok==1)&&(err_oxy<=error_thr) 
    to_end=1; 
  end
   
  
  % draw the results every Ndraw and at the end of simulation 
  if (mod(iter,Ndraw)==0)||(to_end==1)
      
    % add average oxygen and maxmal force to the history vector  
    Nav=Nav+1;  
    avg_oxy(Nav)=mean(mean(oxy));
    max_for(Nav)=max(max(max(forceStr)),max(max(forceTum)));
    
    
    % this is a color map like in teh experiments
    oxy_map=[ 0,255,255;  0,255,255;  0,215,255;  0,174,255;  0,134,255;...
           0, 94,255;  0, 40,255;  0, 13,255;  5, 11,226; 12, 27,183;...
          47, 29, 93;102,  0,  0;168,  0,  0;233,  0,  0;255, 26,  0;...
         255, 64,  0;255,102,  0;255,153,  0;255,178,  0;255,217,  0;...
         255,255,  0;255,255,225]/255;

     
     
    % clear the figure 
    clf
           
    subplot(2,5,[1,2,6,7])  %4 windows for cells and oxygen
      axis([xmin,xmax,ymin,ymax])
      axis equal
      hold on
  
      imagesc(xxg,yyg,oxy(2:Ngx+1,2:Ngy+1)') 
      colorbar
      colormap(oxy_map)
      caxis([0,oxy_max])
      
      viscircles(xyVess,.5*diamVess*ones(Nvess,1),    'color',[0.5,0,0]);
      viscircles(xyStrom,.5*diamStroma*ones(Nstroma,1),'color',[1,1,1]*0.25);
      viscircles(xyTum, .5*diamCell*ones(Ntum,1),     'color', 'k');
   
      if (iter>0)
        quiver(xyTum(:,1),xyTum(:,2),forceTum(:,1),forceTum(:,2),0,'m')
        quiver(xyStrom(:,1),xyStrom(:,2),forceStr(:,1),forceStr(:,2),0,'m')
      end
      
      axis([xmin,xmax,ymin,ymax])
      axis square
      
      title(['[Nv,Nt,Ns]=[',num2str(Nvess),...
          ',',num2str(Ntum),',',num2str(Nstroma),'];   Cell frac=',...
          num2str((Tfrac+Sfrac)*100),'%;   Vasc frac=',num2str(Vfrac*100),'%'])
               

      
      
    subplot(2,5,[3,4,8,9])  %4 windows for oxygen contour
      axis([xmin,xmax,ymin,ymax])
      axis equal
      hold on
  
      contour(xxg,yyg,oxy(2:Ngx+1,2:Ngy+1)','Showtext','on','Linewidth',2) 
      colorbar
      colormap(jet)
      colormap(oxy_map)
      grid minor
   
      axis([xmin,xmax,ymin,ymax])
      axis square
      
      title(['Average O2: ',num2str(mean(mean(oxy))),';   iter=',num2str(iter)])
      
      
      
      
    subplot(2,5,[5,10])  % 2 windows for average oxygen history
       plot(0:Ndraw:(Nav-1)*Ndraw,avg_oxy(1:Nav), 'bo')
       axis([0,Nav*Ndraw,0,60])
       grid on 
       axis square
       
       title(['02 err=',num2str(err_oxy),';  max force: ',num2str(max_for(Nav))])   


    pause(0.1) 
    
    if (toSave==1)&&(mod(iter,Nsave)==0)
      print('-djpeg',[pathFigs,'/fig_V',num2str(round(Vfrac*100)),'_T',...
      num2str(round(Tfrac*100)),'_S',num2str(round(Sfrac*100)),'_',...
      num2str(iter)])
    end
    
  end  % end drawing   
    
    
  % save data
  if (toSave==1)&&(mod(iter,Nsave)==0)
    save([pathData,'/oxy_',num2str(iter),'.txt'],'oxy','-ascii')
    save([pathData,'/xyVess_',num2str(iter),'.txt'],'xyVess','-ascii')
    save([pathData,'/xyTum_',num2str(iter),'.txt'],'xyTum','-ascii')
    save([pathData,'/xyStrom_',num2str(iter),'.txt'],'xyStrom','-ascii')
    avg_oxy_tosave=avg_oxy';
    save([pathData,'/avg_oxy_',num2str(iter),'.txt'],'avg_oxy_tosave','-ascii')
  
    parameters=[seed,xmin,xmax,ymin,ymax,eps,hg,dt,Nav,Nvess,Ntum,...
       Nstroma,Ngx,Ngy,stiff,nu,force_thr,error_thr,Vfrac,Tfrac,Sfrac,...
       diamVess,diamCell,diamStroma,Doxy,oxy_max,oxy_upTum,oxy_upStrom]';
    save([pathDir,'/parameters.txt'],'parameters','-ascii')
  end
end  % end of while loop
  
aver=avg_oxy(Nav);

              
if (toSavelast==1)  
  print('-djpeg',[pathFigs,'/fig_V',num2str(round(Vfrac*100)),'_T',...
        num2str(round(Tfrac*100)),'_S',num2str(round(Sfrac*100))])
   
  save([pathData,'/oxy.txt'],'oxy','-ascii')
  save([pathData,'/xyVess.txt'],'xyVess','-ascii')
  save([pathData,'/xyTum.txt'],'xyTum','-ascii')  
  save([pathData,'/xyStrom.txt'],'xyStrom','-ascii')
  avg_oxy_tosave=avg_oxy';
    
  save([pathData,'/avg_oxy.txt'],'avg_oxy_tosave','-ascii')
    
  parameters=[seed,xmin,xmax,ymin,ymax,eps,hg,dt,Nav,Nvess,Ntum,...
     Nstroma,Ngx,Ngy,stiff,nu,force_thr,error_thr,Vfrac,Tfrac,Sfrac,...
     diamVess,diamCell,diamStroma,Doxy,oxy_max,oxy_upTum,oxy_upStrom]';
  save([pathDir,'/parameters.txt'],'parameters','-ascii')
end
  
  
end % if no_go==1

%----------------------------------------------------------
%----------------------------------------------------------




