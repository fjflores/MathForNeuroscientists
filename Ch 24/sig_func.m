function y = sig_func( p, x )
%compute a sigmoid of the input data x, parameters of the sigmoid are given
%in p
%   y = a + b/(1+exp(-(x-c)/d))

y = p(1) + p(2)./(1 + exp(-(x-p(3))/p(4)));

end

