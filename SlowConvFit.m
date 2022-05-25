tic
h=6.63e-34;
hbar=1.055e-34;
m=39.96*(1.67e-27); % argon mass
L0=0.1559; % Grating-Detector distance
t=160e-9; % Grating substrate thickness
beta=7*(pi/180); % Slit edge angle
Lsd=0.1618; % 2nd Slit-Detector distance
Lss=0.0305; % 1stslit-2nd slit distance
ds2=30e-6; % 2nd slit width
% X1=-0.0053184; % x minimum from signal
% X2=0.0052956; % x maximum from signal
Length=916; % signal pixels
d=257e-9; % grating centre-center separation
s=90e-9; % grating slit width
Ns=20; % number of illuminated slits
A=0.00028; % area of atomic beam signal profile
w1=0.00019; % width of atomic beam signal profile
Ab=0.00046; % area of atomic beam background profile
w1b=0.0053; % width of atomic beam background profile


% Import experimental data
%     xData=importdata('C:\Users\MOT lab\Documents\Data\yData.csv')';
    Data12a=transpose(importdata('C:\Users\MOT lab\Documents\Data\12ms-1.csv'));
    Data12=Data12a(2,:);
    
    % % Use this for the full data set
    data = Data12;
%     x = xData;

    % % Use this to Sample Data by half 
    % data = Data12(2:2:end);
    % x=xData(2:2:end);

    % % Use this to Sample Data by a quarter, 10th etc
    % TrimyD=xData(4:912);
    % TrimD=Data12(4:912);
    % data= TrimD(10:10:end);
    % x=TrimxD(10:10:end);

    % % use this to zoom in on central half of data
%     TrimxD=xData(256:664);
%     TrimD=Data12(256:664);
%     data= TrimD;
%     x=TrimxD;


% xData=importdata('C:\Users\MOT lab\Documents\Data\yData.csv')';
% x = xData;
x=linspace(-0.005327,0.005287,916);
xconv=linspace(-0.0053,0.0053,915);

Res=length(x);

Zeta=linspace(35*(s/2/Res),s/2,Res)';
div=linspace(-0.5,0.5,Res);
VcArray=@(v) arrayfun(@(div0) v(1)+v(2)*div0,div'); % v(1)= average V, v(2) = velocity spread dV 

% Calculate transmission function and slit function

slit=@(y1,V,C3) trapz(Zeta,cos((2*pi*m*V./h).*(y1./L0).*((s./2)-Zeta)).*exp(((j*C3*(10^(-49))*t*cos(beta))./(V.*(Zeta.^3)*hbar)).*((1+(t./2*Zeta)*tan(beta))./(1+(t./Zeta)*tan(beta)).^2)),1);
trans=@(y1,V,C3) (((2./sqrt(h./(m.*V))).*cos((y1)./L0).*slit(y1,V,C3)).*conj((2./sqrt(h./(m.*V))).*cos((y1)./L0).*slit(y1,V,C3))).*...
                (sin((Ns./2)*(2*pi*m*V/h).*(y1)./L0*d)./(sin((1./2)*(2*pi*m*V./h).*(y1)./L0*d))).^2;       
 
% Convolve with detector resolution, numerically integrate variables, normalise
profile=@(x1) ((Ab./(w1b*sqrt(pi./2)))*exp(-2*((x1./w1b).^2))+(A./(w1*sqrt(pi./2)))*exp(-2*((x1./w1).^2)));
convInt=@(x1,v) conv2(trans(x1,VcArray(v),v(3)),profile(xconv),'same');
Int=@(x1,v) trapz(VcArray(v),convInt(x1,v))/v(2); 
finalInt =@(v,x1) (Int(x1,v)/(max(Int(x1,v)))); % normalise the fitted signal

% % % Fitting
vInitial = [13.9,2,8];
vLowerB = [13.60,0.5,4];
vUpperB = [14.20,4,15];
options = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
[v,resnorm,residual,exitflag,output,lamda,J]=lsqcurvefit(finalInt,vInitial,x,data,vLowerB,vUpperB,options);
v
ci = nlparci(v,residual,'jacobian',J)
% fprintf(['The ''trust-region-reflective'' algorithm took %d function evaluations,\n',...
%    'and the ''levenberg-marquardt'' algorithm took %d function evaluations.\n'],...
%    output.funcCount,output2.funcCount)


% % Visualise the sum of least squares for all variables
%  points=8;
%  XV = linspace(13,14.6,points);
%  YD = linspace(1.4,2.2,points);
%  %ZC = linspace(13,17,points)
% for i=17
% for m = 1:points
%     for n = 1:points
%           test(m,n)=sum((data-finalInt([XV(m),YD(n),i],x)).^2);
%       end
% end
% C={'red','green','blue','cyan','magenta','yellow','black','white'};
% hold on
% figure
% % surf(XV,YD,test,'FaceColor',C{i-11})
% surf(XV,YD,test)
% % plot(YV,test)
% end

% Save and plot figures
% v=[13.8827    1.8    12]; % Input paramters for average velocity V, velocity spred dV and C3
fit=[x;finalInt(v,x)]';
fig=figure;
hold on
plot(x,finalInt(v,x))
plot(x,data,'g')
scatter(x,data,5,'*')
print(fig,'C:\Users\MOT lab\Desktop\slow','-dpng');
writematrix(fit,'C:\Users\MOT lab\Desktop\slow.csv');
% 
% halfplot=finalInt(v,x);
% [pks,locs,widths,proms]=findpeaks(halfplot(450:end),x(450:end));
% hold on
% plot(x(450:end),halfplot(450:end))
% plot(locs,pks,'+')
% hold off
toc