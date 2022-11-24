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
a=0.001;%1/sqrt(pi);
t0=5;
b=[0 1];
D0=0.001;
w=2;
g=1;
clrs={'k' ':b' '-.r'};
t=0:0.01:20;
x1=zeros(1,length(t));
x2=zeros(1,length(t));
velocity=zeros(1,length(t));
Ilength=zeros(1,length(t));
    


for j=1:length(b)
            figure
            set(gcf,'color','w');
%set(gcf, 'Position',  [100, 100, 1000, 400])
            for r=1:length(t)
                 z = mu(g,w,t(r),x0,v0)+chol(nearestSPD(sigmaDv2plot(g,w,t(r),t0,a,b(j),d11,d12,d22,D0)))*rand(2,1);
                 x1(r)=z(1);
                 x2(r)=z(2);
            end                            
            plot(real(x1),real(x2),'.b','LineWidth',1)
            xlabel('$x_1$','Interpreter','Latex','FontSize', 10)
            ylabel('$x_2$','Interpreter','Latex','FontSize', 10)
            hold on
            aux=mu(g,w,t,x0,v0);
            plot(real(aux(1,:)),real(aux(2,:)),'k','LineWidth',1)
            %str = sprintf('\\gamma=%.2f', g);
            %sgtitle(str)
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
end
