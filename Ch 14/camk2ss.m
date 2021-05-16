%
% Gabbiani & Cox, Mathematics for Neuroscientists
%
%  camk2ss.m
%
% Compute the bistability curve of Zhabotinsky
%
function camk2ss
close all

pcam = struct('KM',0.4,'KH1',4,'KH2',0.7,...
              'k1',5e-1,'k2',2,'k3',1,'k4',1e-3,...
              'P0',20,'D0',0.05,'I0',.1,...
              'nuCaN',1,'nuPKA',1);

w = [1 16/9 7/3 284/105 25/9 284/105 7/3 16/9 1];

wp(1) = 1;
for j=2:9
    wp(j) = w(j)*wp(j-1);
end
wp = [1 wp];

for Ca=0.05:0.0001:1,
    
    tmp = (Ca/pcam.KH2)^3;
    I = pcam.nuPKA*pcam.I0*(1+tmp)/tmp/pcam.nuCaN;
    
    D = pcam.k4*pcam.D0/(pcam.k3*I + pcam.k4);
    
    tmp = (Ca/pcam.KH1)^4;
    F = tmp/(1+tmp);
    f1 = 10*pcam.k1*F^2;
    f2 = pcam.k1*F;
    f2w = f1*(f2.^[0:9]).*wp./factorial(0:9);    %c_j
    capr = [0 pcam.P0*[1:10].*f2w 0];
    
    capl1 = [-pcam.KM  pcam.k2*D];
    capl2 = [1 f2w];
    capl = conv(capl1,capl2);
    cap = capl - capr;
    
    z = roots(cap);
    rind = find(real(z)>0 & abs(imag(z))<1e-16);
    
    if length(rind)>0
        b = real(z(rind));
        for k=1:length(b)
            b1 = b(k);
            PA = pcam.k2*D/b1 - pcam.KM;
            plot(Ca,PA,'.')
            hold on
        end
    end
    
end
hold off
xlabel('[Ca^{2+}]  (\muM)','fontsize',14)
ylabel('[P]_A   (\muM)','fontsize',14)

return
