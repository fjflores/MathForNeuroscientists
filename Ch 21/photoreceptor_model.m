%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%
%  requires fly_phot9.m

t_v = [1:1500];
rint_ts = zeros(1,1500);
rint_ts(1:500) = 500;
rint_ts(1001:1500) = 500;

c = [-0.8 -0.4 -0.2 -0.1 0.1 0.2 0.4 0.8 1.6 3.2 6.4];

h_f1 = figure;
h_a1 = axes;
for i = 1:length(c)
    rint_ts(501:1000) = 500 + c(i)*500;
    [x_a x_b x_c x_d x_e x_f x_g x_h] = fly_phot9(rint_ts);

    %plot each of the raw variables
    line('Parent',h_a1,'XData',t_v,'YData',x_h,'Color','k');
end;

xlabel(h_a1,'time (ms)');
ylabel(h_a1,'response (mV)');

%print(handles.figure1,'-depsc2','photoreceptor_model.eps');

