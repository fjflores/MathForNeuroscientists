%
%  Gabbiani & Cox, Mathematics for Neuroscientists, 2nd ed
%

h_f = figure; 
h_a1 = subplot(2,2,1);
h_a2 = subplot(2,2,3);
h_a3 = subplot(2,2,2);
h_a4 = subplot(2,2,4);

%uncorrelated
x = -3:0.2:3;
fx = normpdf(x,1,1);
z = fx'*fx;

mesh(h_a1,x,x,z,'EdgeColor','black');
contour(h_a3,x,x,z,'EdgeColor','black');

%correlated
x = -3:0.2:3;
n_x = length(x);
gz = [1 0.8;0.8 1];
gzm = inv(gz);
dgz = det(gz);
c_z = 1/(2*pi*sqrt(dgz));

z = zeros(n_x,n_x);
for i = 1:n_x
    for j = 1:n_x
        v = [x(i) x(j)];
        vz = v - [1 1];
        z(i,j) = c_z*exp(-0.5*vz*gzm*vz');
    end;
end;

mesh(h_a2,x,x,z,'EdgeColor','black');
contour(h_a4,x,x,z,'EdgeColor','black');

