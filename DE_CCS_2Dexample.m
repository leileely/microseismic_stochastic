%% cross-correlation stacking with DE

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
tind = size(recW5,1)*(0:(size(recW5,2)-1));
nnt=size(recW5,1);
recm1z=zeros(nn,mm);

tic;
%% DE parameters
F=0.5;
CR=0.9;
ng=15;%number of generation
np=30;%number of people
%% generate the initial source locations
x=zeros(ng,np);
z=zeros(ng,np);
x(1,:)=round(1+rand(1,np).*(nn-1));
z(1,:)=round(1+rand(1,np).*(mm-1));

ggmax=zeros(1,ng);
gmax=zeros(1,np);



    figure
    subplot(211)
    plot(round(x(1,:)),round(z(1,:)),'k.','markersize',15);
    xlabel('X (m)','fontsize',16);
    ylabel('Z (m)','fontsize',16);
    axis([1 901 1 301]);box on;
    title('iteration no. = 1','fontweight','normal');
    set(gca,'fontsize',16,'xtick',[1:300:901],'xticklabel',10*[],'ytick',[1:100:301],'yticklabel',10*[],'ydir','reverse');


for i=2:ng            
        for p=1:np            
        rr=randperm(np);
        jjj=randperm(2);
        n=rr(1);l=rr(2);m=rr(3);
        if n==p
            n=rr(4);
        elseif l==p
            l=rr(4);
        elseif m==p
            m=rr(4);
        end
        
        if rand<=CR||jjj(1)==1
            x(i,p)=x(i-1,n)+round(F*(x(i-1,l)-x(i-1,m)));
        else
            x(i,p)=x(i-1,p);
        end
        if rand<=CR||jjj(1)==3
            z(i,p)=z(i-1,n)+round(F*(z(i-1,l)-z(i-1,m)));
         else
            z(i,p)=z(i-1,p);
        end
        
         if x(i,p)<1||x(i,p)>nn
            x(i,p)=round(1+rand(1,1).*(nn-1));
        end
               
         if z(i,p)<1||z(i,p)>mm
            z(i,p)=round(1+rand(1,1).*(mm-1));
         end
        
       recm1z(x(i,p),z(i,p))=0;%initialization
        %% CCS stacking process
        ntp = meshgrid(TTPS(:,x(i,p),z(i,p)));
        ntpp = ntp'-ntp;
        ntpp = round(ntpp/dt)+k0; 
        ntpp = min(max(ntpp,1),nnt);
        ntpp = reshape(ntpp,1,s0*s0);
        recm1z(x(i,p),z(i,p)) = sum(recW5(ntpp+tind));
        
        if recm1z(x(i,p),z(i,p))>gmax(p)
            gmax(p)=recm1z(x(i,p),z(i,p));
            xgmax=x(i,p);
            zgmax=z(i,p);
       else
           x(i,p)=x(i-1,p);
           z(i,p)=z(i-1,p);
       end
            
       end
       
       ggmax(i)=max(gmax);


if i==15
    subplot(212)
    plot(round(x(i,:)),round(z(i,:)),'k.','markersize',15);
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

