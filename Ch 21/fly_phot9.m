%
%  Gabbiani & Cox, Mathematics for Neuroscientists
%

function [x_a x_b x_c x_d x_e x_f x_g x_h] = fly_phot9(rint)
%This function is a port of van Hateren's photoreceptor model from Fortran
%to Matlab (directly obtained from the author). In addition, the code 
%has been partially rewritten and commented to make it more intelligible.
     
  ncts=300000; % !max. number of luminance values with 1 sample per ms (5 min)
  ncts_fa=25000; %25 s at a sampling rate of 1 sample per ms
  nt=4096;
  nt2=2*nt;
    
  %%%%%%%%%%%%%%%%%
  %check input data
  %%%%%%%%%%%%%%%%%
  nlen = length(rint);
  if ( nlen == ncts )
    info_str = sprintf('maximal number of values reached %i',ncts);
  else
    info_str = sprintf('read %i values',nlen);
  end;
  disp(info_str);
  
  aver=mean(rint);

  info_str = sprintf('average = %.2f (model assumes average ~600)',aver);
  disp(info_str);

  %%%%%%%%%%%%%%%%%%%%%
  %
  %Model initialization
  %
  %%%%%%%%%%%%%%%%%%%%%
  %final output
  x_h = zeros(1,nlen);

  g_fwd0 = wiener_filter(nt);
  
  %the first filter is 3rd order lowpass, n1 = 3 in vH's original file
  n1 = 3;
  tau1=1.69; %time constant in ms
  [one_m_alpha alpha] = lp(tau1);
    
  %the second filter is 1st order lowpass, n2 = 1 in vH's original file
  n2 = 1;
  tau2=71.8; %time constant in ms
  [one_m_beta beta] = lp(tau2);
  
  %coefficients of the non linearity associated with
  %the 3rd low pass filter
  rk1=0.689;
  rk2=9.07;
  
  %power exponent of 3rd low pass filter
  gcoef=-0.5;	
  %with ncts_fa = 25000, ntaua is equal to 13
  %this is used to determine the number of time constants
  %used to cover the 25 sec of memory using base 4 logarithmic
  %units. Three additional time constants are added on each 
  %side and thus the trailing factor 6.
  ntaua=fix(log10(ncts_fa)/log10(4)+0.5)+6;

  %initialize variables for 3rd filter
  tau_faca = zeros(1,ntaua);	

  ndif=12;
  difa = zeros(1,ndif);

  [o_m_alpha_ms alpha_ms tau_faca ratioa difa] = kern_frac_adap(gcoef,ntaua,ncts_fa,ndif);
  
  %convert to column vector
  o_m_alpha_ms = o_m_alpha_ms(:);
  alpha_ms = alpha_ms(:);
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%
  % initialize and run model
  %%%%%%%%%%%%%%%%%%%%%%%%%%			 
  
  %input steady state value
  prev0=rint(1);
  
  %steady state value after the first divisive
  %non-linearity
  sqrt_prev0=sqrt(prev0);
  
  %steady state value after the second non-linearity
  x_d_st=rtbis(@weber_root,1e-4,1e4,1e-4,rk1,rk2,sqrt_prev0);
  x_f_st=exp(rk2*x_d_st);
  
  %work our way to the steady state-value of the single LP filters
  respdtot=sum(difa(1:ndif))*x_d_st;
  resptot=(x_d_st-respdtot)/ratioa;
  
  tau_factot=sum(tau_faca(1:ntaua));
  
  %steady state value for each of the single LP filters
  rsp_est=resptot/tau_factot;

  %initialize p2mem
  p2mem(1:ndif) = x_d_st;
  
  %save intermediate simulation data
  
  %output of LP1
  x_a = zeros(n1+1,nlen);
  x_a(1:n1+1,1) = prev0;
  
  %output of 1st feedback loop
  x_b = zeros(1,nlen);
  x_b(1) = sqrt_prev0;
  
  %output of LP2
  x_c = zeros(1,nlen);
  x_c(1) = sqrt_prev0;
  
  %output of 2nd feedback loop
  x_d = zeros(1,nlen);
  x_d(1) = x_d_st;
  
  %output of multi-scale LP3 filter
  x_e = zeros(1,nlen);
  x_e(1) = x_d_st;
  
  %output of single low-pass filters
  x_es = zeros(ntaua,nlen);
  x_es(1:ntaua,1) = rsp_est;
  
  %weighted sum ouptput of single low-pass filters
  x_e1 = zeros(1,nlen);
  
  %difference term
  x_e2 = zeros(1,nlen);
  
  %output of 1st static non-linearity
  x_f = zeros(1,nlen);
  x_f(1) = x_f_st;
  
  %output of 2nd static non-linearity
  x_g = zeros(1,nlen);
  x_g(1) = x_d_st./(x_d_st+1);
  
  for i=2:nlen
    
    %feed in the next input
    x_a(1,i) = rint(i);
     
    %compute iteratively the 3rd order LP output
    for j = 2:n1+1
      x_a(j,i) = one_m_alpha*x_a(j-1,i) + alpha*x_a(j,i-1);
    end;
    
    x_c(i) = one_m_beta*x_b(i-1) + beta*x_c(i-1);
    
    %implement the first divisive non-linearity
    x_b(i) = x_a(n1+1,i-1)/x_c(i);
    
    %compute the individual first order LP dynamics
    x_es(1:ntaua,i) = o_m_alpha_ms * x_d(i-1) + alpha_ms .* x_es(1:ntaua,i-1);
    %multi-scale response
    x_e1(i) = ratioa*(tau_faca*x_es(1:ntaua,i));
    
    %correction term, scale the value of x_d ndif terms in the 
    %past by difa and sum the result
    x_e2(i) = difa*p2mem';
    
    %output of the third filter
    x_e(i) = x_e1(i) + x_e2(i);
    
    arg = rk2*x_e(i);
    if (arg > 18) then
      x_f(i) = 1e9;
    else
      %non-linearity following the third filter  
      x_f(i) = exp(rk2*x_e(i));
    end;
    
    %implement the divisive non-linearity following
    %the third order filter after scaling by rk1
    x_d(i) = x_b(i-1)/(rk1*x_f(i));
    
    %final Naka-Rushton non-linearity of the model
    x_g(i) = x_d(i)/(x_d(i)+1);
    
    %save input values for next iteration
    p2mem(ndif:-1:2) = p2mem(ndif-1:-1:1);
    p2mem(1) = x_d(i);
  end;

   
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %
  % convolve with the Wiener filter
  %  
  
  %complex variables used for fft and convolution with
  %wiener filter
  car = zeros(1,nt2);
  car_cts = zeros(1,nt2);
  
  %real variables used for fft and convolution with 
  %wiener filter
  far_cts =  zeros(1,nt);
  far_ctsp = ones(1,nt)*x_g(1);

  jar=1;
  jfar=1;
  
  done = 0; 
  while ( ~done ) 
    if ( jar <= nlen )
      far_cts(jfar)=x_g(jar);
    else
      far_cts(jfar) = 0;
    end;
    
    jfar=jfar+1;
    if (jfar > nt)
      jfar=1;
      
      car_cts(1:nt)=far_ctsp;
      far_ctsp=far_cts;
      
      car(1:nt)=g_fwd0;
      
      car_cts(nt+1:nt2)=far_cts(1:nt);
      car(nt+1:nt2)=0;
      
      %fft of Num Rec in C
      %i.e., nt2*ifft(car_cts) is equivalent to four1(car_cts,nt2,1)
      car_cts = ifft(car_cts);
      car = nt2*ifft(car);
      car=car/nt;

      car=car.*car_cts;
      
      %this is equivalent to four1(car,nt2,-1)
      car = fft(car);
      
      for i=1:nt
        iel=jar+i-nt;
        if (iel <= nlen)
          x_h(iel)=real(car(i+nt));
        end;
      end;
    
      if (jar >= nlen)
        %finished   
        done = 1;
      end;
    end;
    jar=jar+1;
  end

  

%%%%%%%%%%%%%%
%
% Wiener filter
%
%%%%%%%%%%%%%%
function w_f = wiener_filter(nt)

  %sets up the Wiener filter according to vHS, fig. 10, the filter is
  %of the form A(t/tau)^n exp(-t/tau)
  
  w_f = zeros(1,nt);
  
  rsh=0;
  tau=0.535; %time constant, in ms
  n=11; %exponent, n 
  ampl=0.01284; %conversion factor to mV 
  
  for i=1:nt
    arg0=(i+rsh)/tau;
    if (arg0 > 40)
      w_f(i)=0;
    else
      w_f(i)=ampl*(arg0^n)*exp(-arg0);
    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%first order lowpass filter
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [a_0 a_1] = lp(tau)
  %first order low pass filter
  a_1=exp(-1/tau);
  a_0=1-a_1;
  		
function [ar0a ar1a tau_faca ratioa difa] = kern_frac_adap(rcoef,ntaua,ncts_fa,ndif)
  %ndif = 12 and ntaua = 13 are typical 

  %approximates an asymptotic -0.5 power law decay, as in eq. 116
  %of Kasdin
  asy_ma = zeros(1,ncts_fa);
  asy_ma(1)=1;
  for k=2:ncts_fa
    asy_ma(k)=(k-2-rcoef)*asy_ma(k-1)/(k-1);
  end;
  
  %time constants and corresponding weighting factors
  rtau = zeros(1,ntaua);
  dtau = zeros(1,ntaua);

  itau_m3 = (1:ntaua) - 3;
  %time constants on a base 4 logarithmic scale with 3 on each side of
  %the base interval
  rtau=4.^(itau_m3); 
  dtau=4.^(itau_m3); %sampling step in von Schweidler's formula
  rcp_tau=1./rtau;
  %weight factor implementing -1/2 tail according to von Schweidler
  tau_faca=rcp_tau.^(1+rcoef).*dtau;
  
  %corresponding first order exponential LP filters
  ar1a=exp(-rcp_tau);
  ar0a=1-ar1a;

  %first order lp response
  resp_fo   = zeros(1,ntaua);  
  %resp_1 is used to save previous step response; initialized to zero
  resp_1 = zeros(1,ntaua);

  %compound multiscale response to a pulse
  resp_p_ms  = zeros(1,ncts_fa);
  
  %delta function input pulse, first value 1, remaining values zero
  delta_puls = zeros(1,ncts_fa);
  delta_puls(1)=1; 
    
  %compute impulse response of the multi-scale filter
  for i=1:ncts_fa
    %single exponential dynamics
    resp_fo=ar0a*delta_puls(i)+ar1a.*resp_1;
    %multi-scale response
    resp_p_ms(i)=resp_p_ms(i)+resp_fo*tau_faca';
    %save input data for next iteration step
    resp_1=resp_fo;
  end;
  
  %normalizes to approximate as well as possible the Kasdin filter
  %and compute the difference
  nts=50;
  ratioa_range = asy_ma(nts:ncts_fa/2)./resp_p_ms(nts:ncts_fa/2);
  ratioa= mean(ratioa_range);
  difa = asy_ma - ratioa*resp_p_ms;
  
  %compute step response
  step_puls(1:ncts_fa) = ones(1,ncts_fa);
  
  %compound multiscale response to a step
  resp_s_ms  = zeros(1,ncts_fa);
  
  %difference of step from target asymptotic response 
  resp_s_diff = zeros(1,ncts_fa);
  
  %compute step response
  prev0kl = 1;
  resp_1(1:ntaua) = ones(1,ntaua);
  for i=1:ncts_fa
    %single exponential dynamics
    resp_fo=ar0a*step_puls(i)+ar1a.*resp_1;
    %multi-scale response; this is constant 
    resp_s_ms(i)=resp_s_ms(i)+resp_fo*tau_faca';
    %save input data for next iteration step
    resp_1=resp_fo;
    for idif=1:ndif
      iel=i+1-idif;
      if (iel >= 1)
        %response difference; constant
        resp_s_diff(i)=resp_s_diff(i)+difa(idif)*step_puls(iel);
      else
        resp_s_diff(i)=resp_s_diff(i)+difa(idif)*prev0kl;
      end
    end;
  end;
 
  %response, corrected for difference up to ndif time points
  resp_s_ms_corr=ratioa*resp_s_ms+resp_s_diff;
  aver=mean(resp_s_ms_corr);
  %normalize ratioa and difa so that the mean response to a step is one
  ratioa=ratioa/aver;
  difa=difa/aver;
  
  %keep only the first ndif numbers
  difa = difa(1:ndif);
		
function w_r = weber_root(x,rk1,rk2,prev01)

  w_r=rk2*x+log(rk1*x)-log(prev01);

  
function rt_bis = rtbis(func,x1,x2,xacc,rk1,rk2,prev01)
%finds a root of the function func, between x1 and x2 with
%an accuracy xacc

  JMAX=40;
  fmid=func(x2,rk1,rk2,prev01);
  f=func(x1,rk1,rk2,prev01);
  
  if(f*fmid >= 0)
       disp('warning: root must be bracketed in rtbis');
  end

  if(f < 0 )
    rt_bis=x1;
    dx=x2-x1;
  else
    rt_bis=x2;
    dx=x1-x2;
  end

  for j=1:JMAX
    dx=dx*0.5;
    xmid=rt_bis+dx;
    fmid=func(xmid,rk1,rk2,prev01);
    if(fmid <= 0 )
      rt_bis=xmid;
    end;
    if( (abs(dx) < xacc) || (fmid == 0) ) 
      return;
    end;
  end;
  disp('too many bisections in rtbis');
