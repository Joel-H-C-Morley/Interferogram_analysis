tic
% The vectors to create ranges for each variable of each integral: dy, Zeta and V
h=6.626e-34;
hbar=1.055e-34;
m=39.96*(1.67e-27); % argon mass
L0=0.1559; % Grating-Detector distance +/- 0.001
t=160e-9; % Grating substrate thickness +/- 5e-9
beta=7*(pi/180); % Slit edge angle  +/- 0.5
% Lsd=0.1618; % 2nd Slit-Detector distance
% Lss=0.0305; % 1stslit-2nd slit distance
% ds2=30e-6; % 2nd slit width
X1=-0.0053184; % x minimum from signal
X2=0.0052956; % x maximum from signal
% Length=916; % signal pixels
d=257e-9; % grating centre-center separation +/- 5e-9
s=90e-9; % grating slit width +/- 5e-9
Ns=30; % number of illuminated slits
A=0.00028; % area of atomic beam signal profile
w1=0.00020; % width of atomic beam signal profile
Ab=0.00001; % area of atomic beam background profile
w1b=0.0053; % width of atomic beam background profile

% xData=importdata('C:\Users\MOT lab\Documents\Data\yData.csv')';
% x = xData;
trimx=linspace(-1e-3,1e-3,40);
x=linspace(-5.3e-3,5.3e-3,916);

Res=length(x);

% Zeta=linspace(1e-50,s/2,Res)';
Zeta=linspace(40*(s/2/Res),s/2,Res)';
div=linspace(-0.5,0.5,Res);
VcArray=@(v) arrayfun(@(div0) v(1)+v(2)*div0,div'); % v(1)= average V, v(2) = velocity spread dV 

% Calculate transmission function and slit function

slit=@(y1,V,C3) trapz(Zeta,cos((2*pi*m*V./h).*(y1./L0).*((s./2)-Zeta)).*exp(((j*C3*(10^(-49))*t*cos(beta))./(V.*(Zeta.^3)*hbar)).*((1+(t./2*Zeta)*tan(beta))./(1+(t./Zeta)*tan(beta)).^2)),1);
trans=@(y1,V,C3) (((2./sqrt(h./(m.*V))).*cos((y1)./L0).*slit(y1,V,C3)).*conj((2./sqrt(h./(m.*V))).*cos((y1)./L0).*slit(y1,V,C3))).*...
                (sin((Ns./2)*(2*pi*m*V/h).*(y1)./L0*d)./(sin((1./2)*(2*pi*m*V./h).*(y1)./L0*d))).^2;       
      
% slit=@(y1,V,C3) trapz(Zeta,cos((2*pi*m*V./h).*(y1./L0).*((s./2)-Zeta)).*exp(((j*C3*(10^(-49))*t*cos(beta))./(V.*(Zeta.^3)*hbar)).*((1+(t./2*Zeta)*tan(beta))./(1+(t./Zeta)*tan(beta)).^2)),1);
% trans=@(y1,V,C3) (slit(y1,V,C3)).*conj(slit(y1,V,C3));       
            
% Numerically integrate variables, convolve with detector resolution, normalise
% Int=@(x1,v) trapz(VcArray(v),trans(x1,VcArray(v),v(3))); 
%Int2=@(v,x0) arrayfun(@(x1) Int(x1,v)/v(2),x0); % Apply transmission function for each y position
% convInt=@(v,x0) conv(Int2(v,x0),((0.000014/(0.000050*sqrt(pi/2)))*exp(-2*((x/0.000050).^2))),'same');
% convInt=@(x0,v) conv(Int(x0,v),((Ab./(w1*sqrt(pi./2)))*exp(-2*((x./w1b).^2))+(A./(w1*sqrt(pi./2)))*exp(-2*((x./w1).^2))),'same');
% finalInt =@(v,x0) (convInt(x0,v)/(max(convInt(x0,v)))); % normalise the fitted signal
% finalInt =@(v,x0) (Int2(v,x0)/(max(Int2(v,x0)))); % normalise the fitted signal
convInt=@(x1,v) conv2(trans(x1,VcArray(v),v(3)),((Ab./(w1*sqrt(pi./2)))*exp(-2*((x./w1b).^2))+(A./(w1*sqrt(pi./2)))*exp(-2*((x./w1).^2))),'same');
Int=@(x1,v) trapz(VcArray(v),convInt(x1,v))/v(2); 
finalInt =@(x0,v) (Int(x0,v)/(max(Int(x0,v)))); % normalise the fitted signal

% for dv=1:4
% v=[12,dv, 12];
% % output=finalInt(v,x);
% hold on
% plotty=finalInt(x,v);
% % plot(plotty(1,:))
% plot(x(499:end),plotty(499:end))
% ans2(:,dv)=plotty(499:end);
% end
% ans3=ans2/max(max(ans2));
% 
%     Data12a=transpose(importdata('C:\Users\MOT lab\Documents\Data\12ms-1.csv'));
%      Data12=Data12a(2,:);
%     hold on
%     scatter(x,Data12,'+')

% [pks,locs]=findpeaks(Int2(v,x),x,'MinPeakWidth',0.0001);
% plot(locs,pks,'+')
% locs(8)

% for vc =1:1:3
% for dv=1:1:5
%     halfplot=finalInt(x,[13+(1/1),0.0001+(dv*0.4),12]);
%    [pks,locs,widths,proms]=findpeaks(halfplot(450:end),x(450:end));
% %     output=pos;
% pos0(dv)=pks(4+vc);
% pos(dv)=pks(4+vc)-pos0(1);
% end
% hold on
%     plot(pos);
%     output(:,vc)=pos;
% end
% 
% for pk=1:5
% for dv =1:1:51
%     halfplot=finalInt(x,[14,0.0001+(25*0.08),0.645+(dv-1)*0.245]);
%    [pks,locs,widths,proms]=findpeaks(halfplot(499:end),x(499:end));
% %     output=pos;
% pos(dv)=pks(pk);
% pos1(dv)=proms(pk);
% end
% hold on
%    plot(pos);
% %     plot(pos1);
%   output(:,pk)=pos;
%   output1(:,pk)=pos1;
% end

% xtest=linspace(4,10,1000);
% hold on
% plot(xtest,trans(0.002,14,xtest))

hold on
halfplot=finalInt(x,[14,1,12.903]);
[pks,locs,widths,proms]=findpeaks(halfplot(450:end),x(450:end));
plot(x(450:end),halfplot(450:end))
plot(locs,pks,'+');
output=halfplot(450:end);
xout=x(450:end);

% xpp=linspace(0,2.5,150)';
% y0=[-35.75201,-19.67039,-15.26556,-15.2];
% xc=[157.89702,120.29357,131.5167,451.49221];
% w=[105.81181,81.18685,90.00022,82.64639];
% A=[36.1418,20.00004,15.62667,15.5132];
% 
% yp=@(xpp) y0+A.*sin(pi.*(xpp-xc)./w)
% plot(xpp,yp(xpp))
% ans=yp(xpp);

% for i=1:120
% % slope(i)=mean([output(i,7)-output(i,6),output(i,6)-output(i,5),output(i,5)-output(i,4),output(i,4)-output(i,3)]);
% slope(i,1)=  output(i,7)-output(i,6)
% slope(i,2)=  output(i,6)-output(i,5);
% slope(i,3)=  output(i,5)-output(i,4);
% slope(i,4)=  output(i,4)-output(i,3);
% end
% plot(slope)
toc