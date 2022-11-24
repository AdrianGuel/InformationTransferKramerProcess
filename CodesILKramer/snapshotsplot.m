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
c=0.2;
t0=4;
t20=4;
b=0;
d=0;
D0=0.001;
tf1=5;
tf2=5;
w=2;
g=[1];
clrs={'k' ':b' '-.r'};

x1 = -1:0.01:0.5;
x2 = -0.5:0.01:1.2;
[X1,X2] = meshgrid(x1,x2);
X = [X1(:) X2(:)];

for j=1:length(g)
fig = figure;
left_color = [0 0 0];
right_color = [0 0 1];
set(fig,'defaultAxesColorOrder',[left_color; right_color]);
set(fig,'color','w');
set(fig, 'Position',  [100, 100, 800, 400])
    for i=1:length(a)
            t=0:0.01:tf1;
            Ilength=zeros(1,length(t));
            subplot(1,2,1)
            for r=1:length(t)
                if (mod(r,30)==0 && r<1000) || r==1
                    mu=means(g(j),w,t(r),t20,c,d,x0,v0);
                    y = mvnpdf(X,mu',nearestSPD(sigmaDv2plot(g(j),w,t(r),t0,a,b,d11,d12,d22,D0)));
                    y = reshape(y,length(x2),length(x1));
                    contour(real(x1),real(x2),real(y),1);
                    hold on
                    if (mod(r,50)==0 && r<1000) || r==1
                    text(mu(1),mu(2),sprintf('t=%.1f',t(r)))
                    end
                end
                if mod(r,1000)==0
                    y = mvnpdf(X,means(g(j),w,t(r),t20,c,d,x0,v0)',nearestSPD(sigmaDv2plot(g(j),w,t(r),t0,a,b,d11,d12,d22,D0)));
                    y = reshape(y,length(x2),length(x1));
                    contour(real(x1),real(x2),real(y))
                    hold on
                end
            end 
            %axis([-1.5 1.5 -1.5 1.5])
            grid on
            xlabel('$x_1$','Interpreter','Latex','FontSize', 15)
            ylabel('$x_2$','Interpreter','Latex','FontSize', 15)
            hold on
            aux=means(g(j),w,t,t20,c,d,x0,v0);
            plot(real(aux(1,:)),real(aux(2,:)),'k:','LineWidth',1)
            %axis([-01 0.5 -0.5 1.2])
            axis('equal')
            leg1 = legend('$\mathbf{x}\sim N(\langle\mathbf{x}\rangle,\Sigma)$');
            set(leg1,'Interpreter','latex');  
            subplot(1,2,2)
            velocity=EpsilonsDv(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0);
            velocity(isnan(velocity))=0;
            velocitym=Epsilonsmarginals(g(j),w,a(i),b,c,d,t,t0,t20,x0,v0,d11,d12,d22,D0);
            velocitym(isnan(velocitym))=0;
            yyaxis left
            plot(t,real(velocity-(velocitym(1,:)+velocitym(2,:))),clrs{i},'LineWidth',1)
            grid on
            xlabel('$t$','Interpreter','Latex','FontSize', 15)
            ylabel('$\mathcal{E}(t)-(\mathcal{E}_{x_1}(t)+\mathcal{E}_{x_2}(t))$','Interpreter','Latex','FontSize', 15)
            yyaxis right
             plot(t,D(D0,a(i),b,t0,t),'b:','LineWidth',1)
            xlabel('$t$','Interpreter','Latex','FontSize', 15)
            ylabel('$D(t)$','Interpreter','Latex','FontSize', 15)
            leg1 = legend('$\mathcal{E}(t)-\mathcal{E}_m(t)$','$D(t)$');
            set(leg1,'Interpreter','latex');
            axis('equal')
             
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
