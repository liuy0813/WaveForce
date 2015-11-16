% For Tsunami res plot
% 1 plot grid depth of each level
% 2 Arrive plot
close all;clear all;clc
Nestlevel=4;
dt=2;   %seconds
stp=2;  %To save memory
filepath='f:\work\大小鱼山\海啸\CaseTest\';
if exist([filepath,'res'],'dir')
    rmdir([filepath,'res'],'s');
    mkdir([filepath,'res']);
else
    mkdir([filepath,'res']);
end

for ii=1:1 %Nestlevel
    str_id=sprintf('%02i',ii);
    disp('========================')
    disp('Plot Depth Of Each Level')
    disp('========================')
    layer = load([filepath,'layer',str_id,'.dat']);
    layer_x = load([filepath,'layer',str_id,'_x.dat']);
    layer_y = load([filepath,'layer',str_id,'_y.dat']);
    [x,y] = meshgrid(layer_x,layer_y);
    lon0=x(1:stp:end,1:stp:end); lat0=y(1:stp:end,1:stp:end);
    clear x y 
    nx = length(layer_x);ny = length(layer_y);
    clear layer_x layer_x
    depth = reshape(layer,nx,ny);
    layer0=depth(1:stp:end,1:stp:end);
    depth=[];clear depth
    layer=[];clear layer
    layer0(layer0<0)=NaN;
    pcolor(lon0,lat0,layer0');shading flat ;colorbar
    title(['Depth of Level ',str_id])
    figurename = [filepath,'res\','Depth_of_Level_',str_id];
    print(gcf,'-dpng', figurename);
    close;
    disp('========================')
    disp('Plot Zmax Of Level 1 ')
    disp('========================')
    if ii==1
        for jj=3:3
            num0=sprintf('%04i',jj);
            filename=[filepath,'zmax_layer',str_id,'_',num0,'hrs.dat'];
            fid = fopen(filename);
            a = fscanf(fid,'%g',inf); % write all data into a column of matrix a.
            fclose(fid);
            Maxz0 = reshape(a,nx,ny);
            Maxz=Maxz0(1:stp:end,1:stp:end);
            Maxz0=[];clear Maxz0
            m_proj('mercator',...
              'lon',[min(lon0(:)) max(lon0(:))],...
              'lat',[min(lat0(:)) max(lat0(:))]);
            m_pcolor(lon0,lat0,Maxz');shading flat ;
            %hh=contourcmap([0:(4+4)/20:4],'jet','colorbar','on','location','vertical')
            Colortmp=blaq./255;Colortmp(1,:)=[1 1 1];
            colormap(Colortmp);
            colorbar;
            set(gca,'Clim',[0 4]);
            title(['Zmax of Level ',str_id])
            
            Maxz2=Maxz';Maxz_2=Maxz2;
            Maxz2(Maxz_2>0.03)=0;
            Maxz2(Maxz_2<=0.03)=1;
            [numx0,numy0]=size(Maxz2);
            xx0=zeros(1000000,1);yy0=xx0;
            indtmp=1;
            for kk=2:numx0-1
                for hh=2:numy0-1
                    if Maxz2(kk,hh)==0 && layer0(hh,kk)>0
                        if  Maxz2(kk+1,hh)==1 || Maxz2(kk,hh+1)==1 ...
                          ||Maxz2(kk-1,hh)==1 || Maxz2(kk,hh-1)==1 
%                             xx0(indtmp)=lon0(kk,hh);
%                             yy0(indtmp)=lat0(kk,hh);
%                             indtmp=indtmp+1;
%                             Maxz2(kk,hh)=2;
                        end
                    end
                end
            end
            hold on
            m_contour(lon0,lat0,Maxz2,[1 1],'k','Linewidth',1.)
            m_gshhs_l('patch',[0.6 0.6 0.6]);
            m_grid('box','fancy',...
                    'xtick',5,'ytick',5,'tickdir','in',...
                    'fontsize',10);
            figurename = [filepath,'res\','Zmax_of_Level_',str_id,'_',num2str(num0),'Hours'];
            print(gcf,'-dpng', figurename);
%             close;
%             yy0(xx0==0)=[];xx0(xx0==0)=[];
%             xx0=lon0(Maxz2>0.1);
%             yy0=lat0(Maxz2>0.1);
%             xx0(yy0<0)=[];yy0(yy0<0)=[];
%             yy0(xx0<110)=[];xx0(xx0<110)=[];
%             plot(xx0,yy0,'r.')
%             dt = delaunayTriangulation(xx0,yy0);
%             k = convexHull(dt);
%             plot(dt.Points(:,1),dt.Points(:,2), '.', 'markersize',10); hold on;
%             plot(dt.Points(k,1),dt.Points(k,2), 'r'); hold off;

        end   
    end
end