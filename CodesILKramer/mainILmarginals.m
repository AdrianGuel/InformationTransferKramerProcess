%Information length for harmonically bound particle
%Adrian J Guel C
%March 2021

clear all;
close all;
clc

%mu = @(g,w,t,x0,v0) [exp(1).^((-1/2).*g.*t).*(x0.*cosh((1/2).*t.*(g.^2+(-4).*w.^2).^(1/2))+( ...
%  g.^2+(-4).*w.^2).^(-1/2).*(2.*v0+g.*x0).*sinh((1/2).*t.*(g.^2+(-4).* ...
%  w.^2).^(1/2)));exp(1).^((-1/2).*g.*t).*(v0.*cosh((1/2).*t.*(g.^2+(-4).* ...
%  w.^2).^(1/2))+(-1).*(g.^2+(-4).*w.^2).^(-1/2).*(g.*v0+2.*w.^2.*x0).* ...
%  sinh((1/2).*t.*(g.^2+(-4).*w.^2).^(1/2)))];

mucd= @(g,w,t,x0,v0) [exp(1).^((-1).*t.*w).*(x0+t.*(v0+(g+(-1).*w).*x0));exp(1).^((-1).*t.*w) ...
  .*(v0+(-1).*t.*v0.*w+(-1).*t.*w.^2.*x0)];

D= @(D0,a,b,t0,t) D0+b*(exp(-((t-t0).^2)/(a.^2))/(sqrt(pi)*abs(a)));

x0=-0.5;
v0=0.7;
d11=.01;
d12=0;
d21=0;
d22=.01;
a=0.1;%1/sqrt(pi);
c=0.1;
t0=5;
t20=1;
b=0;
d=0;
D0=0.001;
w=1;
g=[1e-6];
clrs={'k' ':b' '-.r'};
t=0:0.01:200;
x1=zeros(1,length(t));
x2=zeros(1,length(t));
velocity=zeros(1,length(t));
velocitym=zeros(1,length(t));
Ilength=zeros(1,length(t));

for j=1:length(g)
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
set(fig,'color','w');
set(fig, 'Position',  [100, 100, 1000, 400])
    for i=1:length(a)
%             for r=1:10:length(t)
%                  z = means(g(j),w,t(r),t20,c,d,x0,v0)+chol(nearestSPD(sigmaDv2plot(g(j),w,t(r),t0,a,b,d11,d12,d22,D0)))*rand(2,1);
%                  x1(r)=z(1);
%                  x2(r)=z(2);
%             end
            velocity=EpsilonsDv(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0);
            velocity(isnan(velocity))=0;
            velocitym=Epsilonsmarginals(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0);
            velocitym(isnan(velocitym))=0;
%             for k=2:length(t)
%                 Ilength(k)=trapz(t(1:k),sqrt(velocity(1:k)));
%             end            
           
        subplot(4,2,1)
         yyaxis left
            plot(t,real(velocity),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$\mathcal{E}(t)$','Interpreter','Latex','FontSize', 10)
            yyaxis right
             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$D(t)$','Interpreter','Latex','FontSize', 10)
         subplot(4,2,3)
          yyaxis left
            plot(t,real(velocitym(1,:)),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$\mathcal{E}_x(t)$','Interpreter','Latex','FontSize', 10)
            yyaxis right
             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$D(t)$','Interpreter','Latex','FontSize', 10)
         subplot(4,2,5)
          yyaxis left
            plot(t,real(velocitym(2,:)),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$\mathcal{E}_v(t)$','Interpreter','Latex','FontSize', 10)
            yyaxis right
             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$D(t)$','Interpreter','Latex','FontSize', 10)
         subplot(4,2,7)
          yyaxis left
            plot(t,real(velocitym(1,:)+velocitym(2,:)-velocity),clrs{i},'LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$\mathcal{E}_x(t)+\mathcal{E}_v(t)-\mathcal{E}(t)$','Interpreter','Latex','FontSize', 10)
            yyaxis right
             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 10)
            ylabel('$D(t)$','Interpreter','Latex','FontSize', 10)             
        subplot(4,2,[2 4 6 8])        
            %plot(real(x1),real(x2),'.b','LineWidth',1)
            %hold on
            %aux=mu(g(j),w,t,x0,v0);
            plot(real(velocitym(1,:)),real(velocitym(2,:)),'k','LineWidth',1)  
            xlabel('$\mathcal{E}_x(t)$','Interpreter','Latex','FontSize', 10)
            ylabel('$\mathcal{E}_v(t)$','Interpreter','Latex','FontSize', 10)            
            % leg1 = legend('$\mathbf{x}\sim N(\langle\mathbf{x}\rangle,\Sigma)$','$\langle\mathbf{x}\rangle$');
            %set(leg1,'Interpreter','latex');             
     end
%     if b==1
%         leg = string(a);
%         for l = 1:length(a)
%            leg(l) = 'a = ' + leg(l);
%         end
%         legend(leg)
%     else
%         leg = 'D_0 = ' + string(D0);
%         legend(leg)
%     end
    str = sprintf('\\gamma=%.2f', g(j));
    sgtitle(str)
end