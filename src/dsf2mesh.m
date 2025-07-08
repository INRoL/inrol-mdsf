function dsf2mesh(objout_path, obj_list, q_list, surf_dir, n1, n2)
% q_list = (Nc, 7)

vi_grid = reshape(1:n1*n2, n1, n2);
fi_1 = vi_grid(2:end, 1:end-1);
fi_2 = vi_grid(1:end-1, 1:end-1);
fi_3 = vi_grid(1:end-1, 2:end);
fi_4 = vi_grid(2:end, 2:end);

fi_u = cat(3, fi_1, fi_2, fi_3);
fi_l = cat(3, fi_3, fi_4, fi_1);
fi = [reshape(fi_u, [], 3); reshape(fi_l, [], 3)];
vnum = size(surf_dir, 2);

mesh_v = [];
mesh_f = [];

for i = 1:length(obj_list)

    [unused, s_res] = obj_list(i).support(q_list(i, :), surf_dir, false, false, false, false);
    
    s = s_res.s;

    mesh_v = [mesh_v; s];
    mesh_f = [mesh_f; fi + vnum*(i-1)];
end

objout = fopen(objout_path, "w");
fprintf(objout, "v %.4f %.4f %.4f\n", mesh_v');
fprintf(objout, "f %d %d %d\n", mesh_f');
fclose(objout);

end