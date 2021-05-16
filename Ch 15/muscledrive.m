%
%  Gabbiani and Cox, Mathematics for Neuroscientists
%
% muscledrive.m
%
% make figures for text from muscle.m
%
% usage:    muscledrive
%
close all
stim = struct('kp',10,'caec',0,'ip3ec',0);
muscle(.005,200,stim)
figure(1)
text(20,1,'(A)','fontsize',20)
figure(2)
text(20,-10,'(B)','fontsize',20)

stim = struct('kp',0,'caec',2,'ip3ec',0);
muscle(.005,200,stim)
figure(3)
text(20,1,'(C)','fontsize',20)
figure(4)
text(20,-15,'(D)','fontsize',20)
