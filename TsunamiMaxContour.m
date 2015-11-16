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
    disp('========================')
    disp('Plot Zmax Of Level 1 ')
    disp('========================')            
    filename0=[filepath,'zmax_layer01.dat'];
    fid = fopen(filename0);
    a = fscanf(fid,'%g',inf); % write all data into a column of matrix a.
    fclose(fid);
    Maxzlayer01 = reshape(a,nx,ny);
    Maxz01=Maxzlayer01(1:stp:end,1:stp:end);
    Maxz0=[];clear Maxzlayer01
    m_proj('mercator',...
           'lon',[min(lon0(:)) max(lon0(:))],...
           'lat',[min(lat0(:)) max(lat0(:))]);
    m_pcolor(lon0,lat0,Maxz01');shading flat ;
    %hh=contourcmap([0:(4+4)/20:4],'jet','colorbar','on','location','vertical')
    Colortmp=blaq./255;Colortmp(1,:)=[1 1 1];
    colormap(Colortmp);
    colorbar;
    set(gca,'Clim',[0 4]);
    title(['Zmax of Level ',str_id])
    for jj=1:15
        num0=sprintf('%04i',jj);
        filename=[filepath,'zmax_layer',str_id,'_',num0,'hrs.dat'];
        fid = fopen(filename);
        a = fscanf(fid,'%g',inf); % write all data into a column of matrix a.
        fclose(fid);
        Maxz0 = reshape(a,nx,ny);
        Maxz=Maxz0(1:stp:end,1:stp:end);
        Maxz0=[];clear Maxz0
        Maxz2=Maxz';Maxz_2=Maxz2;
        Maxz2(Maxz_2>0.03)=0;
        Maxz2(Maxz_2<=0.03)=1;
        hold on
        m_contour(lon0,lat0,Maxz2,[1 1],'k','Linewidth',1.)
     end
     m_gshhs_l('patch',[0.6 0.6 0.6]);
     m_grid('box','fancy',...
            'xtick',5,'ytick',5,'tickdir','in',...
            'fontsize',10);
     figurename = [filepath,'res\','Zmax_of_Level_',str_id,'_',num2str(num0),'Hours'];
%             print(gcf,'-dpng', figurename);
end