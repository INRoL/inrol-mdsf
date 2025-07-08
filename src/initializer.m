function para = initializer(V, maxp)

n_th1 = 31;
n_th2 = 2*(n_th1 - 1) + 1;
th1_range = pi / 2;
th2_range = pi;
th1 = linspace(-th1_range, th1_range, n_th1);
th2 = linspace(-th2_range, th2_range, n_th2);
cos_th1 = cos(th1);
sin_th1 = sin(th1);
cos_th2 = cos(th2);
sin_th2 = sin(th2);

directions = [];
for i = 2:(n_th1-1)
  for j = 2:(n_th2)
    d_x = cos_th1(i)*cos_th2(j);
    d_y = cos_th1(i)*sin_th2(j);
    d_z = sin_th1(i);
    d = reshape([d_x, d_y, d_z], [1, 3]);
    directions = [directions; d];
  end
end
directions = [directions; 0, 0, 1; 0, 0, -1];
rot_directions = directions * rotx(13) * roty(55) * rotz(71);





n_th1 = 101;
n_th2 = 2*(n_th1 - 1) + 1;

th1 = linspace(-pi/2,pi/2,n_th1);
th2 = linspace(0,2*pi,n_th2);
[VV,UU] = meshgrid(th1,th2);
cosV = cos(VV); % error needs to be clamped at k*pi/2
sinV = sin(VV);
cosU = cos(UU);
sinU = sin(UU); % error needs to be clamped at k*pi
cosV(:,[1 n_th1]) = 0;
sinU([1 n_th2],:) = 0;

x_n = cosV.*cosU;
y_n = cosV.*sinU;
z_n = sinV;
n1 = size(x_n, 1);
n2 = size(x_n, 2);
surfdir = reshape(cat(3, x_n, y_n, z_n), [], 3);


% if k = 1/V, interp = 1/V
alpha = zeros(1, V);

phi = (1 + sqrt(5))/2;
dodeca = [1.0, 1.0, 1.0; 1.0, 1.0, -1.0; 1.0, -1.0, 1.0; 1.0, -1.0, -1.0; -1.0, 1.0, 1.0; -1.0, 1.0, -1.0; -1.0, -1.0, 1.0; -1.0, -1.0, -1.0; 0.0, phi, 1/phi; 0.0, phi, -1/phi; 0.0, -phi, 1/phi; 0.0, -phi, -1/phi; 1/phi, 0.0, phi; 1/phi, 0.0, -phi; -1/phi, 0.0, phi; -1/phi, 0.0, -phi; phi, 1/phi, 0.0; phi, -1/phi, 0.0; -phi, 1/phi, 0.0; -phi, -1/phi, 0.0];
init_vert = (dodeca * 0.5 / phi)';




%% Parameters
%----- Constants
para.V = V;
para.maxp = maxp;

%----- Initial shape
para.init_v = init_vert; % (3,V), width=1
if V > 20
    tmp_vert = rand(3, V-20);
    tmp_vert = 2*tmp_vert - 1;
    para.init_v = [para.init_v, tmp_vert];
end
if V < 20
    para.init_v = rand(3, V)*2 - 1;
end
para.init_a = alpha; % (1,V)
para.init_p = 0;

%----- Directions
para.directions = directions';
para.rot_directions = rot_directions';

%----- Visualization parameters
para.surf_n1 = n1;
para.surf_n2 = n2;
para.surf_directions = surfdir';
para.verbose = true;
para.showfig = true;

%----- Optimization parameters
para.maxiter = 30;
para.maxfunc = 300;

%----- quaternion
para.quat_id = [0,0,0,1,0,0,0];

end