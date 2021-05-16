%
% Gabbiani and Cox, Mathematics for Neuroscientists
%
% four2b.m  
%
% Illustrate aspects of the 2-dimensional Fourier Transform
% on a brain slice - contained in  mri_image.mat  courtesy of
% Dr. M. S. Beauchamp, Dept. of Neurosurgery, 
% Baylor College of Medicine
%
%  usage:   four2b
%

close all

load('mri_image');

cam = mbmri144;

figure(1)                               % plot original image
x = linspace(-10,10,314);
y = x;
imagesc(x,y,real(cam))
colormap('gray')
axis equal
axis tight
xlabel('y   (cm)','fontsize',14)
ylabel('z   (cm)','fontsize',14)

figure(2)                               % plot its Fourier Transform
fcam = fft2(double(cam));
sfcam = fftshift(fcam);

%plot absolute value of Fourier coefficients 
%and normalize peak coefficient to 1
masfcam = abs(sfcam)/max(abs(sfcam(:)));
kx = linspace(-157/20,157/20,314);
ky = kx;
%the last parameter sets the color scale to show 90 percent of pixels 
imagesc(kx,ky,masfcam,[0.001 0.0029])
colormap(gray)
colorbar
hold on
cutoff = 15;
plot(cutoff*[-1 1 1 -1 -1]/20,cutoff*[-1 -1 1 1 -1]/20,'r','linewidth',1)
hold off
axis equal
axis tight
xlabel('k_y   (1/cm)','fontsize',14)
ylabel('k_z   (1/cm)','fontsize',14)

figure(3)     % plot image recovered from undersampled Fourier Transform
fcamunder = fcam(1:2:end,1:2:end);
camunder = ifft2(fcamunder);
imagesc(x,y,real(camunder))
colormap('gray')
axis equal
axis tight
xlabel('y   (cm)','fontsize',14)
ylabel('z   (cm)','fontsize',14)

figure(4)   % plot image recovered from Low Pass Filtered Fourier Transform
sfcamlo = zeros(314);
range = 157-cutoff:157+cutoff;
sfcamlo(range,range) = sfcam(range,range);
camlo = ifft2(ifftshift(sfcamlo));
imagesc(x,y,real(camlo))
colormap('gray')
axis equal
axis tight
xlabel('y   (cm)','fontsize',14)
ylabel('z   (cm)','fontsize',14)

figure(5)  % plot image recovered from High Pass Filtered Fourier Transform
sfcamhi = sfcam;
sfcamhi(range,range) = zeros(length(range));
camhi = ifft2(ifftshift(sfcamhi));
imagesc(x,y,real(camhi))
colormap('gray')
axis equal
axis tight
xlabel('y   (cm)','fontsize',14)
ylabel('z   (cm)','fontsize',14)

