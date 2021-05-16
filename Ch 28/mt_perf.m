%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%

%uncomment to plot a typical psychometric function
%alpha = 6.1;
%beta = 1.17;

%c = 0.1:0.1:100;
%p = 1 - 0.5*exp( -(c/alpha).^beta );

%figure; semilogx(c,p);


%plot the distributions of spikes in the two conditions
%for 3 coherence levels

c_vals = [3.2 6.4 12.8 25.6 51.2 100];

t_obs = 0.5; %half a second
m_p = (0.265*c_vals + 23.2)*t_obs;
m_a = (-0.072*c_vals + 23.2)*t_obs;

s_p = sqrt(1.5*m_p);
s_a = sqrt(1.5*m_a);

for i = 3:2:5
    figure;
    dx = 6*s_p(i)/100;
    x = m_p(i)-3*s_p(i):dx:m_p(i)+3*s_p(i);
    y = normpdf(x,m_p(i),s_p(i));
    plot(x,y);
    hold on;
 
    dx = 6*s_a(i)/100;
    x = m_a(i)-3*s_a(i):dx:m_a(i)+3*s_a(i);
    y = normpdf(x,m_a(i),s_a(i));
    plot(x,y,'r');
    xlabel('number of spikes');
    ylabel('probability density');

end;

figure; hold on;
pc_vect = zeros(size(c_vals));
pc_2afc_vect = zeros(size(c_vals));
for i = 1:length(c_vals)
    %compute the roc curve
    mu_y = (m_p(i) - m_a(i))/s_a(i);
    s_y = s_p(i)/s_a(i);

    rho_0 = 2*s_y^2*log(s_y);

    l_val_vect = -3:0.05:3;
    pfa_vect = zeros(size(l_val_vect));
    pd_vect = zeros(size(l_val_vect));

    for j = 1:length(l_val_vect)
        l_val = l_val_vect(j);
        rho = rho_0 + 2*s_y^2*l_val;

        d_y2 = s_y^2*(mu_y^2+rho)-rho;
    
        if ( d_y2 > 0 )
            d_y = sqrt(s_y^2*(mu_y^2+rho)-rho);

            y_p = (-mu_y + d_y)/(s_y^2-1);
            y_m = (-mu_y - d_y)/(s_y^2-1);

            pfa_vect(j) = normcdf(y_m,0,1) + (1-normcdf(y_p,0,1));
            pd_vect(j) = normcdf(y_m,mu_y,s_y) + (1-normcdf(y_p,mu_y,s_y));
        else
            pfa_vect(j) = NaN;
            pd_vect(j) = NaN;
        end;
    
    end;

    inds1 = find(~isnan(pfa_vect));
    inds2 = find(~isnan(pd_vect));

    if ( ~isempty(find(inds1 ~= inds2)) )
        disp('error not the same NaN indices');
        return;
    end;

    [pfa_incr, ix] = sort(pfa_vect(inds1));
    pfa_incr = [0 pfa_incr 1];

    pd_incr = pd_vect(inds1(ix));
    pd_incr = [0 pd_incr 1];

    if ( (i == 3) | (i == 5) )
        plot(pfa_incr,pd_incr);
    end;
    
   
    %figure; plot(pfa_incr,pd_incr);
    
    err_incr = 0.5*pfa_incr + 0.5*(1-pd_incr);
    
    %figure; plot(pfa_incr,err_incr);
    
    pc_vect(i) = 1 - min(err_incr);

    %compute the area under the roc curve
    pc_2afc_vect(i) = trapz(pfa_incr,pd_incr);
end;
xlabel('P_F_A');
ylabel('P_D');


figure; semilogx(c_vals,pc_vect);
hold on; semilogx(c_vals,pc_2afc_vect,'r');

alpha = 20;
beta = 1.2;
p_psycho = 1 - 0.5*exp( -(c_vals/alpha).^beta );
semilogx(c_vals,p_psycho,'g');
xlabel('coherence (%)');
ylabel('probability correct');

%s_plot = 10;
%dx = 2*s_plot*s_y/1000;
%x = mu_y-s_plot*s_y:dx:mu_y+s_plot*s_y;
%y = normpdf(x,mu_y,s_y);
%plot(x,y);

%hold on;

%dx = 2*s_plot/1000;
%x = -s_plot:dx:s_plot;
%y = normpdf(x,0,1);
%plot(x,y,'r');