clc
clear 
close all
%%参数设置
z=3300;              % z传输距离,单位m
N=1000;             % N采样点数目
w0=0.012;            % 束腰半径,单位m
a=0.1;               % 光束的指数衰减因子
Cn2=5*10.^-15;       % 湍流结构常数
lambda=1550*10.^-9;  % 波长,单位m
k=2*pi/lambda;

%普通高斯光束参数
z1 =0;
n_i=1.0;                                   % 背景空间折射率
k0=k*n_i;                                  % 背景空间的波数
ZR=k0*w0^2/2;                              % 瑞利距离                                        
w=w0*sqrt(1+(z1/ZR)^2);                    % 传播到z处的束宽
R=ZR*(z1/ZR+ZR/z1);                        % 等相位面曲率半径
Phi=atan(z1/ZR);                           % 相位因子
%普通高斯光束参数

bushu=16;            % 生成15个等间距相位屏
deltz=z/bushu;       % 设置相位屏间距,单位m
M = 1;
m0 =3;               % 涡旋拓扑荷数
d = 0.066;           % 横向位移参数,单位m
%r1=(2.^0.5).*d;      % 主环半径,单位m 
r1= 0.1;
r0 = (0.4229.*(k.^2).*Cn2.*deltz).^(-3/5); %大气相干长度
L0=10;               % 外尺度,单位m
l0=0.01;             % 内尺度,单位m
L = 1;               % 输入屏尺寸,单位m
L2 = 1;              % 相位屏尺寸,单位m
subh = 3;            % 次谐波次数
delta=L/N;           % 步长，即空间采样间隔，单位m
delta1=L2/N;         % 步长，即空间采样间隔，单位m
x=(-N/2:N/2-1)*delta;
y=x;
[X, Y]=meshgrid(x,y);
[theta,r] = cart2pol(X,Y);
%r=(X.^2+Y.^2).^0.5;
%%参数设置
%X = r.*cos(theta);
%Y = r.*sin(theta);

%%高斯光束电场表达式
E0=(1/(1+1i*z1./ZR)).*exp(-(r.^2./w0^2)./(1+1i*z1./ZR)).*exp(1i*k0.*z1);
%%高斯光束电场表达式

figure (1)
I0 =abs(E0).^2;
pcolor(x,y,I0);  
title("源场光场");
colormap("jet")
shading interp;
axis on;
axis square;  %输出正方形图像
xlabel('{\itx}(m)');
ylabel('{\ity}(m)');

% figure (2)
% phase1=(angle(E0)*180/pi);
% pcolor(X,Y,phase1); 
% colormap("Gray");
% shading interp;
% axis on;
% axis square;  %输出正方形图像
% colormap("jet")


%%未加湍流条件下的传输
x1=linspace(-L./2,L./2,N).*N; % *N就给出频域坐标
y1=linspace(-L./2,L./2,N).*N;
[kethi,nenta]=meshgrid(x1,y1);
U0=E0;
H=exp(1i.*k.*z.*(1-(lambda.^2).*(kethi.^2+nenta.^2)./2));  % 传递函数H
G0=fftshift(fft2(U0));                 % 衍射面上光场的傅里叶变换
G=G0.*H;                               % 光场的频谱与传递函数相乘
U=ifft2(G);                            % 在观察屏上的光场分布
I1=U.*conj(U);                         % 在观察屏上的光强分布
figure (3)
pcolor(x,y,I1); 
title("理想传输后的光场");
colormap("jet");
colorbar;
axis on;
shading interp;
axis square;  %输出正方形图像
xlabel('{\itx}(m)');
ylabel('{\ity}(m)');


figure (4)
phase2=(angle(U)*180/pi);
pcolor(X,Y,phase2); 
colormap("Gray");
shading interp;
axis on;
axis square;  %输出正方形图像
colormap("jet")
%%未加湍流条件下的传输

%%加湍流条件下的传输
UT = E0;
del_f=1/(N*delta1);              % 频率间隔
fx=(-N/2:N/2-1)*del_f;           % 横向频率
[kx ,ky]=meshgrid(2*pi*fx);      % 生成网格坐标矩阵
[th ,ka]=cart2pol(kx, ky);       % 直极转化
km=5.92/l0;
k0=2*pi/L0;
PSD_phi=0.033*Cn2*exp(-(ka/km).^2)./(ka.^2+k0.^2).^(11/6);%修正的von karman大气湍流谱模型
PSD_phi(N/2+1,N/2+1)=0;
cn=2.*pi.*(k.^2).*deltz.*PSD_phi.*(2.*pi.*del_f).^2;% 相位扰动的幅度

sum =0;
%for n = 1:2

h=waitbar(0,'计算中，请等待...');
for l=1:bushu

 %% 高频部分
 %生成大小为 N x N 的高斯随机数矩阵，然后将其加上相应幅度的湍流相位扰动
 %最后进行逆傅里叶变换 (ift2) 得到空间域中的相位畸变。
 phz_hi=ift2((randn(N)+1i*randn(N)).*sqrt(cn),1);
 phz_hi=real(phz_hi);
 %% 高频部分

 %% 低频补偿
phz_lo=zeros(size(phz_hi));
for p=1:subh
del_fp=1/(3^p*L2);
fx1=(-1:1)*del_fp;
[kx1, ky1]=meshgrid(2*pi*fx1);
[th1 ,k1]=cart2pol(kx1,ky1);
km=5.92/10;
k0=2*pi/L0;  %outscale frequency
PSD_phi1=0.033*Cn2*exp(-(k1/km).^2)./(k1.^2+k0.^2).^(11/6);
PSD_phi1(2,2)=0;
cn1=2*pi*k.^2*deltz.*PSD_phi1.*(2*pi*del_fp).^2;
cn1=(randn(3)+1i*randn(3)).*sqrt(cn1);
SH=zeros(N);
    for ii=1:9
    SH=SH+cn1(ii)*exp(1i*(kx1(ii)*X+ky1(ii)*Y)) ;
    end
    phz_lo=phz_lo+SH;
 end
  %% 低频补偿

phz_lo=real(phz_lo)-mean(real(phz_lo(:)));
phz=phz_hi+phz_lo;
figure (5)
pcolor(kx,ky,phz); 
colormap("jet")
axis on;
shading interp;
axis square;  %输出正方形图像
colorbar;
title("低频补偿后的大气湍流随机相位屏");
% xlabel('{\itx}(m)');
% ylabel('{\ity}(m)');

x2=linspace(-L/2,L./2,N).*N; % *N就给出频域坐标
y2=linspace(-L/2,L./2,N).*N;
[kethi2,nenta2]=meshgrid(x2,y2);

%%判断是否是第一次循环，第一次循环时不叠加随机相位屏
   if l == 1
       G1=fftshift(fft2(UT));                 %衍射面上光场的傅里叶变换
   else
       G1 = G1.*exp(1i.*phz);
   end
%%判断是否是第一次循环，第一次循环时不叠加随机相位屏

H=exp(1i*k*deltz.*(1-(lambda.^2).*(kethi2.^2+nenta2.^2)./2)); %传递函数H
G1=G1.*H;                                % 光场的频谱与传递函数相乘
Ut=ifft2(G1);                            % 在观察屏上的光场分布
    if l ~= max(bushu)
    I2=Ut.*conj(Ut); 
    figure (6)
    pcolor(x,y,I2);
    nu = l-1;
    title("经"+nu+"次大气湍流传输后的光场");
    colormap("jet")
    axis on;
    shading interp;
    axis square;  %输出正方形图像
    axis image;
    xlabel('{\itx}(m)');
    ylabel('{\ity}(m)');
    end
       shijian=num2str(fix(l/bushu*100));
       waitbar(l/bushu,h,['请等待，已完成',shijian,'%']);
end
close(h);
%%加湍流条件下的传输

%sum = sum + Ut;
%end
%U = sum/2;

%an = 1./((2.*pi).^0.5).*int(Ut.*exp(-1i.*n.*r2),r2,0,2.*pi);

figure (7)
phase3=(angle(Ut)*180/pi);
pcolor(X,Y,phase3); 
colormap("Gray");
shading interp;
axis on;
axis square;  %输出正方形图像
colormap("jet")

I3=Ut.*conj(Ut);                          %在观察屏上的光强分布
figure (8)
pcolor(x,y,I3); 
nu2 = bushu -1;
title("经"+nu2+"次大气湍流传输后的光场");
colormap("jet")
shading interp;
axis square;  %输出正方形图像
axis image;
colorbar;
xlabel('{\itx}(m)');
ylabel('{\ity}(m)');