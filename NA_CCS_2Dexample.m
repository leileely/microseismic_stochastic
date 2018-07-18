%% cross-correlation stacking with NA

clear;clc;

load steiner2drec;
clear recvx;
dt=0.001;
s0=101;
k0=2501;
nn=901;
mm=301;

load TTSTEINERP;%traveltime table
TTPS=permute(TTPS,[3,1,2]);


%% cross-correlation
recW5=zeros(2*(k0-1)+1,s0*s0);
for ii=1:s0
    for jj=ii+1:s0
        recW5(:,ii+(jj-1)*s0)=xcorr(recvz(:,ii),recvz(:,jj));
    end 
end
nnt=size(recW5,1);
recm1z=zeros(1,nn*mm);

tic;
%% NA parameters
np=30;
numb=10;
nump=2;
ng=8;
N=np;
M=np/numb;
%% generate the initial source locations
xmin = diag([1,1])*ones(nump,np); 
xmax = diag([nn,mm])*ones(nump,np); 
x = xmin+(xmax-xmin).*rand(nump,np); 

IMAGE=zeros(1,np);
ggmaxna=zeros(1,ng);
range=[1,nn;1,mm];


X0=round(x);
xx = sub2ind([nn,mm],X0(1,:),X0(2,:));
tind = size(recW5,1)*(0:(size(recW5,2)-1));
 
for p=1:np;
    ntp = meshgrid(TTPS(:,X0(1,p),X0(2,p)));
    ntpp = ntp'-ntp;
    ntpp = round(ntpp/dt)+k0; 
    ntpp = min(max(ntpp,1),nnt);
    ntpp = reshape(ntpp,1,s0*s0);
    recm1z(xx(:,p)) = sum(recW5(ntpp+tind));
    IMAGE(p)=sum(recW5(ntpp+tind));
end
[~,INDEX]=sort(IMAGE,'descend');
ggmaxna(1)=IMAGE(INDEX(1));
[gmax,maxpos]=max(recm1z(xx(:)));
xgmax=X0(:,maxpos)*ones(1,np);


    figure
    subplot(211)
    [vx,vy]=voronoi(X0(1,:),X0(2,:));
    plot(X0(1,:),X0(2,:),'k.','markersize',8);
    hold on;
    plot(vx,vy,'k-','linewidth',1.5)%
    xlabel('X (m)','fontsize',16);
    ylabel('Z (m)','fontsize',16);
    axis([1 901 1 301]);box on;
    title('iteration no. = 1','fontweight','normal');
    set(gca,'fontsize',16,'xtick',[1:300:901],'xticklabel',10*[],'ytick',[1:100:301],'yticklabel',10*[],'ydir','reverse');
%% NA iterations
for i=2:ng    
       
        
        %%updating rules of NA
        x=zeros(nump,np);
        for n=1:numb
        x1=X0(:,INDEX(n));
        x2=X0(:,[1:INDEX(n)-1,INDEX(n)+1:N]);
        x0=zeros(nump,M);
        for m=1:M
            if m==1
                x0(:,m)=x1;
            else
                x0(:,m)=x0(:,m-1);
            end 
            
            Dk1=sum((x0(2,m)-x1(2)).^2);
            Dk2=sum((x0(1,m)-x1(1)).^2);
            Dj1=sum(((repmat(x0(2,m),1,N-1)-x2(2,:)).^2),1);
            Dj2=sum(((repmat(x0(1,m),1,N-1)-x2(1,:)).^2),1);
            Dk=[Dk1;Dk2];
            Dj=[Dj1;Dj2];
            for ii=1:nump
                  
                xji=(repmat(x1(ii),1,N-1)+x2(ii,:)+(Dk(ii)-Dj(ii,:))./(repmat(x1(ii),1,N-1)-x2(ii,:)))/2;
                r1=max([range(ii,1),xji(xji<x0(ii,m))]);
                r2=min([range(ii,2),xji(xji>x0(ii,m))]);
                x0(ii,m)=r1+(r2-r1)*rand(1,1);                
                
            end
        end
        x(:,(n-1)*M+1:n*M)=x0;
        end
        
        X=round(x);
        xx = sub2ind([nn,mm],X(1,:),X(2,:));
        recm1z(xx(:))=0;%initialization
        image=zeros(1,np);
        %% CCS stacking process
        for p=1:np;
        ntp = meshgrid(TTPS(:,X(1,p),X(2,p)));
        ntpp = ntp'-ntp;
        ntpp = round(ntpp/dt)+k0; 
        ntpp = min(max(ntpp,1),nnt);
        ntpp = reshape(ntpp,1,s0*s0);
        recm1z(xx(:,p)) = sum(recW5(ntpp+tind));
        image(p)=sum(recW5(ntpp+tind));
        end
        
        IMAGE=[IMAGE,image];
        X0=[X0,X];
        N=N+np;
        
       
        [~,INDEX]=sort(IMAGE,'descend');
        [gmax,maxpos]=max(recm1z);
        [xg1,xg2]=ind2sub([nn,mm],maxpos);
        xgmax=[xg1;xg2]*ones(1,np);
        ggmaxna(i)=gmax;

if i==8
    subplot(212)
    [vx,vy]=voronoi(X0(1,:),X0(2,:));
    plot(X0(1,:),X0(2,:),'k.','markersize',8);
    hold on;
    plot(vx,vy,'k-','linewidth',1.5)%
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
plot(max(ggmaxna)-ggmaxna);

