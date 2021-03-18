%Information transfer for harmonically bound particle
% when \omega\rightarrow 0 or free brownian motion
%Adrian J Guel C
%January 2021

clear all;
close all;
clc

mu = @(g,w,t,x0,v0) [exp(1).^((-1/2).*g.*t).*(x0.*cosh((1/2).*t.*(g.^2+(-4).*w.^2).^(1/2))+( ...
  g.^2+(-4).*w.^2).^(-1/2).*(2.*v0+g.*x0).*sinh((1/2).*t.*(g.^2+(-4).* ...
  w.^2).^(1/2)));exp(1).^((-1/2).*g.*t).*(v0.*cosh((1/2).*t.*(g.^2+(-4).* ...
  w.^2).^(1/2))+(-1).*(g.^2+(-4).*w.^2).^(-1/2).*(g.*v0+2.*w.^2.*x0).* ...
  sinh((1/2).*t.*(g.^2+(-4).*w.^2).^(1/2)))];

mucd= @(g,w,t,x0,v0) [exp(1).^((-1).*t.*w).*(x0+t.*(v0+(g+(-1).*w).*x0));exp(1).^((-1).*t.*w) ...
  .*(v0+(-1).*t.*v0.*w+(-1).*t.*w.^2.*x0)];

D= @(D0,a,b,t0,t) D0+b*(exp(-((t-t0).^2)/(a.^2))/(sqrt(pi)*abs(a)));

x0=-0.5;
v0=0.7;
d11=.01;
d12=0;
d21=0;
d22=.01;
a=0.1;
t0=1;
b=0;
D0=0.001;
w=1;
g=[1e-6 1 1.9999 3];
clrs={'k' ':b' '-.r'};
t=0:0.001:20;
x1=zeros(1,length(t));
x2=zeros(1,length(t));

for j=1:length(g)
    figure
    set(gcf,'color','w');
    set(gcf, 'Position',  [100, 100, 1000, 400])
    for i=1:length(a)
        if g(j)==2*w
            S=sigmacdDv(g(j),w,t,t0,a(i),b,d11,d12,d22,D0);
            for r=1:length(t)
                 z = mucd(g(j),w,t(r),x0,v0)+chol(nearestSPD(sigmacdDvplot(g(j),w,t(r),t0,a,b,d11,d12,d22,D0)))*rand(2,1);
                 x1(r)=z(1);
                 x2(r)=z(2);        
            end
        else
            S=sigmaDv2(g(j),w,t,t0,a(i),b,d11,d12,d22,D0);
            for r=1:10:length(t)
                 z = mu(g(j),w,t(r),x0,v0)+chol(nearestSPD(sigmaDv2plot(g(j),w,t(r),t0,a,b,d11,d12,d22,D0)))*rand(2,1);
                 x1(r)=z(1);
                 x2(r)=z(2);        
            end            
        end
        
        Hx=0.5*(1+log(2*pi*S(1,:)));
        %dHxdt = gradient(Hx(:)) ./ gradient(t(:));
        dHxdt=0.5*(gradient(S(1,:)) ./ gradient(t(:)'))./S(1,:);    
        Hv=0.5*(1+log(2*pi*S(4,:)));
        dHvdt = gradient(Hv(:)) ./ gradient(t(:));
        Txv=dHvdt'+g(j)-(D(D0,a(i),b,t0,t)-w^2.*S(2,:))./S(4,:);
        
        subplot(3,3,1)
            plot(t,real(Hx),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$H_x(t)$','Interpreter','Latex','FontSize', 10)
            hold on
         subplot(3,3,2)
            plot(t,real(Hv),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$H_v(t)$','Interpreter','Latex','FontSize', 10)
            hold on           
        subplot(3,3,4)
            plot(t,real(dHxdt),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$\frac{dH_x(t)}{dt}$','Interpreter','Latex','FontSize', 10)
            hold on
        subplot(3,3,5)
            plot(t,real(dHvdt),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$\frac{dH_v(t)}{dt}$','Interpreter','Latex','FontSize', 10)
            hold on            
        subplot(3,3,7)
            plot(t,real(dHxdt),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$T_{v\to x}$','Interpreter','Latex','FontSize', 10)
            hold on
        subplot(3,3,8)
            plot(t,real(Txv),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$T_{x\to v}$','Interpreter','Latex','FontSize', 10)
            hold on
        subplot(3,3,[3 6])
            plot(real(x1),real(x2),'.b','LineWidth',1)
            xlabel('$x_1$','Interpreter','Latex','FontSize', 10)
            ylabel('$x_2$','Interpreter','Latex','FontSize', 10)
            hold on
            aux=mu(g(j),w,t,x0,v0);
            plot(real(aux(1,:)),real(aux(2,:)),'k','LineWidth',1)
            hold on
        subplot(3,3,9)
            plot(t,D(D0,a(i),b,t0,t),'b','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$D(t)$','Interpreter','Latex','FontSize', 10)
            hold on                
    end
    if b==1
        leg = string(a);
        for l = 1:length(a)
           leg(l) = 'a = ' + leg(l);
        end
        legend(leg)
    else
        leg = 'D_0 = ' + string(D0);
        legend(leg)
    end
    str = sprintf('\\gamma=%.2f', g(j));
    sgtitle(str)
end