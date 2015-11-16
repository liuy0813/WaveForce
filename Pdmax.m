%
% 波流力计算  2015-11-3 14:46:39
% 根据海港水文规范（2013） 式 8.3.1中的公式
%
close all;clear all;clc
%%
%
% 基础数据
%
d=10.72;  %建筑物前水深
g=9.8;
pi=3.141592653;
T=5.12;         %Period of the wave
H=2.09;         %Height of the wave
L=38.49 ;        %Wave length
Cd=1.2;         %速度力系数,1.2 For Circle and 2.0 for squre
Cm=2.0;         %惯性力系数,2.0 For Circle and 2.2 For squre
Gama=10.25;     %海水容重，γ
D=1.5;          %桩直径
A=pi.*D^2./4;   %桩的截面积
omega=2*pi./T;
yita=0.65*H;    %From Fig.8.3.2-1
aa=0.62;        %波面变形受波高和水深的影响
% dx=[5,5,5,5,5,5,5,5]; %各排桩之间的间距
% dx0=[0,6.6,13.2,20.1,27,33.6,40.2]; %zeros(length(dx)+1,1);
dx0=[0,2.4,4.5,6.6,9,18,27,29.4,31.5,33.6,36];
nt=12;  %将一个波分为12段
outputdir='./FeiYunJiang/';
if exist(outputdir,'dir')
    rmdir(outputdir,'s');
    mkdir(outputdir);
else
    mkdir(outputdir);    
end
wavefile=[outputdir,'WaveEta.txt'];
waveout=[outputdir,'waveout.txt'];
WaveForce=[outputdir,'WaveForce.txt'];
WaveForcePhase=[outputdir,'WaveForcePhase.txt'];
%
%  特定波长和周期下，第一排桩t=0时，各排群桩各处桩位的时刻
%  各排桩距离第一排桩的距离,第一排桩tt=0时，其他各排桩处的时间
%
%  cos(k*x-omega*t) 其中k=2*pi/L x为距离第一排桩的位置 sin 的同理
% 
phase0=zeros(length(dx0),1);
for jj=1:length(dx0)
    phase0(jj)=dx0(jj)*2*pi/L;
end
%
%  一个周期分为nt份，输出一排桩不同时刻处波面高度
%
eta=zeros(length(phase0),nt);
for ii=1:nt
    eta(:,ii)=H.*aa.*cos(phase0-omega.*(ii-1)*2*pi/nt);
end
fid=fopen(wavefile,'w');
fprintf(fid,'%10s','Phase');
for ii=1:length(dx0)
    fprintf(fid,'%10.3f',dx0(ii));
end
fprintf(fid,'\n');
for jj=1:nt
    fprintf(fid,'%10.3f',(jj-1)*2*pi/nt);
    for ii=1:length(dx0)
        fprintf(fid,'%10.3f',eta(ii,jj)+d);
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%
%%
%
%  ==== 8.3.1 作用于底面高度z处柱体全断面上与波向平行的正向力  ====
%  ====                   速度分力+惯性分力
%
% ==== 8.3.2 作用于整个柱体高度上的最大速度分力+最大惯性分力
%    五点法的z,包括z=0 d 和 d+η
%    wt=0°时，η=max(η) ,wt=270°时，η=max(η)-H./2
%    五点法算出来的是五个点上的受力，PDmax算出的是四段上的受力情况
%
z=[0,(d-0)./3,(d-0)*2./3,d,d+yita];
dz=zeros(length(dx0),nt,length(z));
p1=zeros(length(dx0),nt,length(z));PD=zeros(length(dx0),nt,length(z));
%
%  Personal Understand
%  这里算的应该是不同时刻（每个桩都有不同的波面高度，eta），不同的波面高度作用下的波浪力
%  将一个周期分为12个部分，第一根桩有12种不同的起始波面高度
%  这里有问题 PD趋势正确，但p1并不正确,把eta(ii,jj)直接换成H之后 趋势正确!
%  五点法
for jj=1:nt
    for ii=1:length(dx0)
        for kk=1:length(z)
%             z0=z(kk); %(eta(ii,jj)+d)*(kk-1)/4; %z(kk); %eta(ii,jj)+d; %波面高度+水深
            dz(ii,jj,kk)=(eta(ii,jj)+d)*(kk-1)/4; %%  波面高度+水深 为了重现报告中的结果
            z0=dz(ii,jj,kk);
            uu=pi*H*cosh(2*pi*z0/L)*cos(phase0(ii)-omega.*(jj-1)*2*pi/nt)/(sinh(2*pi*d/L)*T);
            pupt=2*(pi^2)*H*cosh(2*pi*z0/L)*sin(phase0(ii)-omega.*(jj-1)*2*pi/nt)/(T^2*sinh(2*pi*d/L));
            p1(ii,jj,kk)=Gama*Cm*A*pupt/g;
            PD(ii,jj,kk)=0.5*Gama*Cd*D*uu*abs(uu)/g;
        end
    end
end
psum0=p1+PD; 
%
% 求总力
%
SumForces=zeros(length(dx0),nt);
for jj=1:nt %  只需要0时刻的值
    for ii=1:length(dx0)
        zz=squeeze(dz(ii,jj,:));
        Forces=squeeze(psum0(ii,jj,:));
        dzz=diff(zz);
        SumForces(ii,jj)=sum((Forces(1:end-1)+Forces(2:end)).*0.5.*dzz);  %梯形面积
    end
end
%
% 输出0时刻五点法各排单桩波浪力计算表  单位kN
%
fid=fopen(WaveForce,'w');
%  第一行
fprintf(fid,'%10s','桩位');
for ii=1:length(dx0)
    fprintf(fid,'%10i',ii);
end
fprintf(fid,'\n');
% 第二行 波面高度
fprintf(fid,'%10s','波面高度');
for ii=1:length(dx0)
    fprintf(fid,'%10.3f',eta(ii,1)+d);
end
fprintf(fid,'\n');
%
% 第三行---第七行 五点法  单桩波浪力
%
for jj=1:5
    fprintf(fid,'%10i   ',jj);
    for ii=1:length(dx0)
        fprintf(fid,'%10.3f',psum0(ii,1,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'%10s ','总力');
for ii=1:length(dx0)
    fprintf(fid,'%10.3f',SumForces(ii,1));
end
fclose(fid);
%
% 输出不同相位单桩波浪力计算表
%
fid=fopen(WaveForcePhase,'w');
fprintf(fid,'%10s','桩位');
for ii=1:length(dx0)
    fprintf(fid,'%10i',ii);
end
fprintf(fid,'\n');
for jj=1:nt
    fprintf(fid,'%10.3f',(jj-1)*2*pi/nt);
    for ii=1:length(dx0)
        fprintf(fid,'%10.3f',SumForces(ii,jj));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%====到此====结束=====
return
%
%  其实我看不懂这一行是要干嘛，只是照着改写的
%
psum1=sum(psum0(:,:,1:end-1)+psum0(:,:,2:end),3).*(eta+repmat(d,[length(dx0) nt]))./8;
fid=fopen(wavefile,'a');
fprintf(fid,'%10s','Phase');
for ii=1:length(dx0)
    fprintf(fid,'%10.3f',dx0(ii));
end
fprintf(fid,'\n');
for jj=1:nt
    fprintf(fid,'%10.3f',(jj-1)*2*pi/nt);
    for ii=1:length(dx0)
        fprintf(fid,'%10.3f',psum1(ii,jj));
    end
    fprintf(fid,'\n');
end
fclose(fid);
%%
%
% PDmax=Cd*Gama*D*H^2*K1/2
% PImax=Cm*Gama*A*H*K2/2
% 五点,共四段
%
PDmax=zeros(length(z)-1,1);PImax=zeros(length(z)-1,1);
for ii=1:length(z)-1
    z1=z(ii);z2=z(ii+1);
    K1=(4*pi*z2/L-4*pi*z1/L+sinh(4*pi*z2/L)-sinh(4*pi*z1/L))./(8*sinh(4*pi*d/L));
    K2=(sinh(2*pi*z2/L)-sinh(2*pi*z1/L))/cosh(2*pi*d/L);
    PDmax(ii)=Cd*Gama*D*H^2*K1/2;
    PImax(ii)=Cm*Gama*A*H*K2/2;
end  
%
%作用于整个柱体高度上任何相位时的正向水平总波浪力
%
P=zeros(length(dx0),nt);
for jj=1:nt
    for ii=1:length(dx0) 
        P(ii,jj)=sum(PDmax).*cos(omega.*(jj-1)*2*pi/nt-phase0(ii)).*abs(cos(omega.*(jj-1)*2*pi/nt-phase0(ii)))-... !!- -->+ because sin(-x)=sin(x)
                  sum(PImax).*sin(phase0(ii)-omega.*(jj-1)*2*pi/nt);
    end
end
%%
return
%%
%
% ==== PDmax 和 PImax对z1断面的力矩MDmax和MImax
%
K3=zeros(length(z)-1,1);K4=zeros(length(z)-1,1);
MDmax=zeros(length(z)-1,1);MImax=zeros(length(z)-1,1);
for ii=1:length(z)-1
    z1=z(ii);z2=z(ii+1);
    K3part01=(pi^2*(z2-z1)^2)/(4*L^2);
    K3part02=pi*(z2-z1)*sinh(4*pi*z2/L)/(8*L);
    K3part03=(cosh(4*pi*z2/L)-cosh(4*pi*z1/L))/32;
    K3part04=sinh(4*pi*d/L);
    K3(ii)=(K3part01+K3part02+K3part03)/K3part04;
    K4part01=2*pi*(z2-z1)/L*sinh(2*pi*z2/L);
    K4part02=cosh(2*pi*z2/L)-cosh(2*pi*z1/L);
    K4part03=cosh(2*pi*d/L);
    K4(ii)=(K4part01-K4part02)/K4part03;
    MDmax(ii)=Cd*Gama*D*H^2*L*K3(ii)/2/pi;
    MImax(ii)=Cm*Gama*A*H*L*K4(ii)/4/pi;
end

