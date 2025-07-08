function obj = mesh2dsf(obj_path, q, para, out_path)
%% Mesh Loading
mesh_data = readObj(obj_path);
vin = mesh_data.v';

maxv = max(vin, [], 2);
minv = min(vin, [], 2);
orig = (maxv + minv)/2; % (3,1)

scale = 1 / max(maxv - minv);
vin = (vin-orig)*scale;

%% Initial parameters
con_nor = para.rot_directions; % (3,D)
con_dist = max(vin'*con_nor, [], 1); % (1,D)

obj_init = DSF3d(para.init_v.*0.1, para.init_a, para.init_p);

%% Optimization
display = "iter";
options = optimoptions("lsqnonlin", "Algorithm", "trust-region-reflective", "Display", display, ...
    "MaxIterations", para.maxiter, "MaxFunctionEvaluations", para.maxfunc, ...
    "SpecifyObjectiveGradient", true);

obj = fitting(obj_init, q, con_nor, con_dist, options, para.verbose, true);

obj = obj.scale(1/scale).translate(orig);

%% Visualization
if para.showfig == true
    faces = mesh_data.f.v;
    vin = vin/scale + orig;

    f = figure();
    f.Position = [0, 0, 1000, 1000];
    trimesh(faces, vin(1, :), vin(2, :), vin(3, :), "facecolor", "g", "facealpha", 0.2, "linestyle", "none");

    hold on;
    xlabel("x");
    ylabel("y");
    zlabel("z");
    xlim([-0.5, 0.5]);
    ylim([-0.5, 0.5]);
    zlim([-0.5, 0.5]);
    pbaspect([1, 1, 1]);
    box on
    grid on;

    color = "b";
    alpha = 0.5;
    edgealpha = 0.2;
    show_v = true;
    show_t = true;

    obj.draw(q, para.surf_directions, para.surf_n1, para.surf_n2, color, alpha, edgealpha, show_v, show_t);

    titstr = strcat("p: ", num2str(obj.p));
    title(titstr);
end

if strlength(out_path) > 1
    dsf2mesh(out_path, obj, q, para.surf_directions, para.surf_n1, para.surf_n2);
end

end







