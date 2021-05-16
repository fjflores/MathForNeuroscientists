%
% LNTcurves.m
% 
% Gabbiani & Cox, Mathematics for Neuroscientists
%
% Ca functionals from
% Jaffe et al, A model for dendritic Ca2+ accumulation in
% hippocampal pyramidal neurons based on fluorescence imaging measurements,
% J Neurophysiol, 71(3):1065-1077, 1994.
%
% usage: LNTcurves(L,N,T)
%
% e.g.
%   L = struct('a_L',15.69,'b_L',81.5,'c_L',0.29,'d_L',10.86)
%   N = struct('a_N',0.19,'b_N',19.88,'c_N',0.046,'d_N',20.73,...
%   'e_N',1.6e-4,'f_N',48.46,'g_N',39)
%   T = struct('a_T',0.2,'b_T',19.26,'c_T',0.009,'d_T',22.03,...
%   'e_T',1e-6,'f_T',16.26,'g_T',29.79)
%


function LNTcurves(L,N,T)

V = -80:.1:40;

alpha_m = L.a_L*(L.b_L-V)./(exp((L.b_L-V)/10)-1);
beta_m = L.c_L*exp(-V/L.d_L);
tau_m = 1./(alpha_m+beta_m);
m_infty = alpha_m.*tau_m;
figure(1)
subplot(5,1,1)
[ax,h1,h2] = plotyy(V,m_infty,V,tau_m);
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','r')
set(h1,'color','k')
set(h2,'color','r')
ylabel(ax(1),'m_{\infty,L}','fontsize',14,'color','k')
ylabel(ax(2),'\tau_{m,L}','fontsize',14,'color','k')
box off

alpha_m = N.a_N*(N.b_N-V)./(exp((N.b_N-V)/10)-1);
beta_m = N.c_N*exp(-V/N.d_N);
tau_m = 1./(alpha_m+beta_m);
m_infty = alpha_m.*tau_m;
alpha_h = N.e_N*exp(-V/N.f_N);
beta_h = 1./(1+exp((N.g_N-V)/10));
tau_h = 1./(alpha_h+beta_h);
h_infty = alpha_h.*tau_h;
subplot(5,1,2)
[ax,h1,h2] = plotyy(V,m_infty,V,tau_m);
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','r')
set(ax(1),'xticklabel',[])
set(ax(2),'xticklabel',[])
set(h1,'color','k')
set(h2,'color','r')
ylabel(ax(1),'m_{\infty,N}','fontsize',14,'color','k')
ylabel(ax(2),'\tau_{m,N}','fontsize',14,'color','k')
box off
subplot(5,1,3)
[ax,h1,h2] = plotyy(V,h_infty,V,tau_h);
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','r')
set(h1,'color','k')
set(h2,'color','r')
ylabel(ax(1),'h_{\infty,N}','fontsize',14,'color','k')
ylabel(ax(2),'\tau_{h,N}','fontsize',14,'color','k')
box off

alpha_m = T.a_T*(T.b_T-V)./(exp((T.b_T-V)/10)-1);
beta_m = T.c_T*exp(-V/T.d_T);
tau_m = 1./(alpha_m+beta_m);
m_infty = alpha_m.*tau_m;
alpha_h = T.e_T*exp(-V/T.f_T);
beta_h = 1./(1+exp((T.g_T-V)/10));
tau_h = 1./(alpha_h+beta_h);
h_infty = alpha_h.*tau_h;
subplot(5,1,4)
[ax,h1,h2] = plotyy(V,m_infty,V,tau_m);
set(ax(1),'xticklabel',[])
set(ax(2),'xticklabel',[])
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','r')
set(h1,'color','k')
set(h2,'color','r')
ylabel(ax(1),'m_{\infty,T}','fontsize',14,'color','k')
ylabel(ax(2),'\tau_{m,T}','fontsize',14,'color','k')
box off
subplot(5,1,5)
[ax,h1,h2] = plotyy(V,h_infty,V,tau_h);
set(ax(1),'ycolor','k')
set(ax(2),'ycolor','r')
set(h1,'color','k')
set(h2,'color','r')
ylabel(ax(1),'h_{\infty,T}','fontsize',14,'color','k')
ylabel(ax(2),'\tau_{h,T}','fontsize',14,'color','k')
xlabel('V   (mV)','fontsize',14,'color','k')
box off


