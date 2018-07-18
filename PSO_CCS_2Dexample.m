%% cross-correlation stacking (CCS) with PSO

clear;clc;

load steiner2drec;
clear recvx;
dt=0.001;
s0=101;
k0=2501;
nn=901;
mm=301;

load TTSTEINERP;%traveltime table
TTPS=reshape(TTPS,nn*mm,s0);

%% cross-correlation
recW5=zeros(2*(k0-1)+1,s0,s0);
for ii=1:s0
    for jj=1:s0
        recW5(:,jj,ii)=xcorr(recvz(:,ii),recvz(:,jj));
    end 
end
clear recvz;

nnt=size(recW5,1);
recm1z=zeros(1,nn*mm);%imaging values

tic;

%% PSO parameters
c1=2;
c2=2;
nump=2;%number of parameters (x,z)
ng=15;%number of generations
np=30;%number of particles
w=1:-(1-0.3)/ng:0.3;%0.7
%% generate the initial source locations
xmin = diag([1,1])*ones(nump,np); 
xmax = diag([nn,mm])*ones(nump,np); 
vmax = diag([nn/3,mm/3])*ones(nump,np); 
x = xmin+(xmax-xmin).*rand(nump,np); 
xx = sub2ind([nn,mm],round(x(1,:)),round(x(2,:)));
v = zeros(nump,np);
xpmax=x;

pmax=zeros(1,np);
ggmax=zeros(1,ng);

%% the first iteration to get the gmax

     for ii=1:s0
         for jj=1:s0                            
                ntp=round((TTPS(xx,ii)-TTPS(xx,jj))/dt)+k0;              
                ntp=min(max(ntp,1),nnt);                 
                pmax=pmax+recW5(ntp,jj,ii)'.^2;                           
                recm1z(xx(:))=recm1z(xx(:))+recW5(ntp,jj,ii)'.^2;
          end
     end

gmax=max(pmax);
ggmax(1)=gmax;
[gmax,maxpos]=max(recm1z(xx(:)));
xgmax=x(:,maxpos)*ones(1,np);

%% GIF-part1
% filename='SteinerPSOexample.gif';
% plot(round(x(1,:)),round(x(2,:)),'k.','markersize',15);
% set(gca,'fontsize',16,'xtick',[1:300:901],'xticklabel',10*[0:300:900],'ytick',[1:100:301],'yticklabel',10*[0:100:300]);
% xlabel('X (m)','fontsize',16);
% ylabel('Z (m)','fontsize',16);
% axis([1 901 1 301]);box on;pause(1);
% frame = getframe(1);        
% im = frame2im(frame);        
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
%%


    figure
    subplot(211)
    plot(round(x(1,:)),round(x(2,:)),'k.','markersize',15);
    xlabel('X (m)','fontsize',16);
    ylabel('Z (m)','fontsize',16);
    axis([1 901 1 301]);box on;
    title('iteration no. = 1','fontweight','normal');
    set(gca,'fontsize',16,'xtick',[1:300:901],'xticklabel',10*[],'ytick',[1:100:301],'yticklabel',10*[],'ydir','reverse');



%% generations/iterations
for i=2:ng   
       
        
        %%updating rules of PSO
        v=w(i-1)*v+c1*rand(nump,np).*(xpmax-x)+c2*rand(nump,np).*(xgmax-x);
        v=min(max(v,-vmax),vmax);
        x=x+v;
        ind = ((x < xmin) + (x > xmax)) == 1;
        xrand = xmin+(xmax-xmin).*rand(nump,np); 
        x(ind) = xrand(ind);   
        
        xx = sub2ind([nn,mm],round(x(1,:)),round(x(2,:)));
        recm1z(xx(:))=0;%initialization
        %% CCS stacking process
        for ii=1:s0            
                 for jj=1:s0            
                ntp=round((TTPS(xx,ii)-TTPS(xx,jj))/dt)+k0;              
                ntp=min(max(ntp,1),nnt);                                      
                recm1z(xx(:))=recm1z(xx(:))+recW5(ntp,jj,ii)'.^2;
                end    
        end
        
        ind = recm1z(xx(:))>pmax;
        pmax(ind) = recm1z(ind);
        xpmax(:,ind)=x(:,ind);
        
        [gmax,maxpos]=max(recm1z);
        [xg1,xg2]=ind2sub([nn,mm],maxpos);
        xgmax=[xg1;xg2]*ones(1,np);
        ggmax(i)=gmax;
       

%% GIF-part2
% plot(round(x(1,:)),round(x(2,:)),'k.','markersize',15);
% set(gca,'fontsize',16,'xtick',[1:300:901],'xticklabel',10*[0:300:900],'ytick',[1:100:301],'yticklabel',10*[0:100:300]);
% xlabel('X (m)','fontsize',16);
% ylabel('Z (m)','fontsize',16);
% axis([1 901 1 301]);box on;pause(1);
% frame = getframe(1);        
% im = frame2im(frame);        
% [imind,cm] = rgb2ind(im,256);
% imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0);%,'LoopCount',inf
%%

if i==15
    subplot(212)
    plot(round(x(1,:)),round(x(2,:)),'k.','markersize',15);
    xlabel('X (m)','fontsize',16);
    ylabel('Z (m)','fontsize',16);
    axis([1 901 1 301]);box on;
    title('iteration no. = 15','fontweight','normal');
    set(gca,'fontsize',16,'xtick',[1:300:901],'xticklabel',10*[],'ytick',[1:100:301],'yticklabel',10*[],'ydir','reverse');
    set(gcf,'pos',[100 100 1000 700]);
end
end

toc;
recm1z=reshape(recm1z,nn,mm);  

ind=find(recm1z==max(recm1z(:)));
[xgmax,zgmax]=ind2sub([nn,mm],ind);

figure%show 2d slice of the result
imagesc((recm1z.^2./max(recm1z(:).^2))')
colormap(jet)
hc=colorbar;
minn=-0.2;
maxx=1;
caxis([minn maxx])
set(hc,'ytick',[minn maxx],'yticklabel',{'low' 'high'},'fontsize',16)
set(hc,'pos',[0.928 0.145 0.02 0.78]);
set(gca,'xtick',[1:100:901],'xticklabel',[0:1000:9000])
set(gca,'ytick',[1:50:301],'yticklabel',[0:500:3000])
set(gca,'fontsize',16);
set(gcf,'pos',[100 100 1000 400])
xlabel('X (m)');
ylabel('Z (m)','fontsize',16);
axis([1 901 1 301])

figure%the converging line
plot(max(ggmax)-ggmax);

