%
% ����������  2015-11-3 14:46:39
% ���ݺ���ˮ�Ĺ淶��2013�� ʽ 8.3.1�еĹ�ʽ
%
close all;clear all;clc
%%
%
% ��������
%
d=10.72;  %������ǰˮ��
g=9.8;
pi=3.141592653;
T=5.12;         %Period of the wave
H=2.09;         %Height of the wave
L=38.49 ;        %Wave length
Cd=1.2;         %�ٶ���ϵ��,1.2 For Circle and 2.0 for squre
Cm=2.0;         %������ϵ��,2.0 For Circle and 2.2 For squre
Gama=10.25;     %��ˮ���أ���
D=1.5;          %׮ֱ��
A=pi.*D^2./4;   %׮�Ľ����
omega=2*pi./T;
yita=0.65*H;    %From Fig.8.3.2-1
aa=0.62;        %��������ܲ��ߺ�ˮ���Ӱ��
% dx=[5,5,5,5,5,5,5,5]; %����׮֮��ļ��
% dx0=[0,6.6,13.2,20.1,27,33.6,40.2]; %zeros(length(dx)+1,1);
dx0=[0,2.4,4.5,6.6,9,18,27,29.4,31.5,33.6,36];
nt=12;  %��һ������Ϊ12��
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
%  �ض������������£���һ��׮t=0ʱ������Ⱥ׮����׮λ��ʱ��
%  ����׮�����һ��׮�ľ���,��һ��׮tt=0ʱ����������׮����ʱ��
%
%  cos(k*x-omega*t) ����k=2*pi/L xΪ�����һ��׮��λ�� sin ��ͬ��
% 
phase0=zeros(length(dx0),1);
for jj=1:length(dx0)
    phase0(jj)=dx0(jj)*2*pi/L;
end
%
%  һ�����ڷ�Ϊnt�ݣ����һ��׮��ͬʱ�̴�����߶�
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
%  ==== 8.3.1 �����ڵ���߶�z������ȫ�������벨��ƽ�е�������  ====
%  ====                   �ٶȷ���+���Է���
%
% ==== 8.3.2 ��������������߶��ϵ�����ٶȷ���+�����Է���
%    ��㷨��z,����z=0 d �� d+��
%    wt=0��ʱ����=max(��) ,wt=270��ʱ����=max(��)-H./2
%    ��㷨���������������ϵ�������PDmax��������Ķ��ϵ��������
%
z=[0,(d-0)./3,(d-0)*2./3,d,d+yita];
dz=zeros(length(dx0),nt,length(z));
p1=zeros(length(dx0),nt,length(z));PD=zeros(length(dx0),nt,length(z));
%
%  Personal Understand
%  �������Ӧ���ǲ�ͬʱ�̣�ÿ��׮���в�ͬ�Ĳ���߶ȣ�eta������ͬ�Ĳ���߶������µĲ�����
%  ��һ�����ڷ�Ϊ12�����֣���һ��׮��12�ֲ�ͬ����ʼ����߶�
%  ���������� PD������ȷ����p1������ȷ,��eta(ii,jj)ֱ�ӻ���H֮�� ������ȷ!
%  ��㷨
for jj=1:nt
    for ii=1:length(dx0)
        for kk=1:length(z)
%             z0=z(kk); %(eta(ii,jj)+d)*(kk-1)/4; %z(kk); %eta(ii,jj)+d; %����߶�+ˮ��
            dz(ii,jj,kk)=(eta(ii,jj)+d)*(kk-1)/4; %%  ����߶�+ˮ�� Ϊ�����ֱ����еĽ��
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
% ������
%
SumForces=zeros(length(dx0),nt);
for jj=1:nt %  ֻ��Ҫ0ʱ�̵�ֵ
    for ii=1:length(dx0)
        zz=squeeze(dz(ii,jj,:));
        Forces=squeeze(psum0(ii,jj,:));
        dzz=diff(zz);
        SumForces(ii,jj)=sum((Forces(1:end-1)+Forces(2:end)).*0.5.*dzz);  %�������
    end
end
%
% ���0ʱ����㷨���ŵ�׮�����������  ��λkN
%
fid=fopen(WaveForce,'w');
%  ��һ��
fprintf(fid,'%10s','׮λ');
for ii=1:length(dx0)
    fprintf(fid,'%10i',ii);
end
fprintf(fid,'\n');
% �ڶ��� ����߶�
fprintf(fid,'%10s','����߶�');
for ii=1:length(dx0)
    fprintf(fid,'%10.3f',eta(ii,1)+d);
end
fprintf(fid,'\n');
%
% ������---������ ��㷨  ��׮������
%
for jj=1:5
    fprintf(fid,'%10i   ',jj);
    for ii=1:length(dx0)
        fprintf(fid,'%10.3f',psum0(ii,1,jj));
    end
    fprintf(fid,'\n');
end
fprintf(fid,'%10s ','����');
for ii=1:length(dx0)
    fprintf(fid,'%10.3f',SumForces(ii,1));
end
fclose(fid);
%
% �����ͬ��λ��׮�����������
%
fid=fopen(WaveForcePhase,'w');
fprintf(fid,'%10s','׮λ');
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
%====����====����=====
return
%
%  ��ʵ�ҿ�������һ����Ҫ���ֻ�����Ÿ�д��
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
% ���,���Ķ�
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
%��������������߶����κ���λʱ������ˮƽ�ܲ�����
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
% ==== PDmax �� PImax��z1���������MDmax��MImax
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

