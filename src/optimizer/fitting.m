function [obj, minc] = fitting(obj_init, q, con_nor, con_dist, options, verbose, search_p)

phi_init = obj_init.get_phi();
objective = @(var)shape_cost(var, q, con_nor, con_dist);
if verbose == true
    tic;
end
phi = lsqnonlin(objective, phi_init', [], [], options);
if verbose == true
    dt = toc;
    disp(['[Shaping] Optimization Time: ', num2str(dt), ' sec']);
end

obj = DSF3d(phi');


%% Search integer p

if (search_p)
    obj_init = obj;
    p = obj_init.p;
    pc_init = get_pcost(obj_init.get_phi(), p, q, con_nor, con_dist);
    
    p_list = [];
    if (floor(p) - 1 >= 2)
        p_list = [p_list, floor(p) - 1];
    end
    if (floor(p) >= 2 && floor(p) ~= p)
        p_list = [p_list, floor(p)];
    end
    if (ceil(p) <= 40 && ceil(p) ~= p)
        p_list = [p_list, ceil(p)];
    end
    if (ceil(p) + 1 <= 40)
        p_list = [p_list, ceil(p) + 1];
    end
    
    c_list = [];
    obj_list = [];
    if (floor(p) == p || ceil(p) == p)
        c_list = [c_list, pc_init];
        obj_list = [obj_list, obj_init];
    end

    for i = 1:length(p_list)
        phi_init = obj_init.get_phi();
        phi_init(end) = [];
        
        objective = @(var)shape_cost_p(var, p_list(i), q, con_nor, con_dist);
        
        if verbose == true
            tic;
        end
        phi = lsqnonlin(objective, phi_init', [], [], options);
        if verbose == true
            dt = toc;
            disp(['[Search p: ', num2str(p_list(i)), '] Optimization Time: ', num2str(dt), ' sec']);
        end
        obj = DSF3d([phi; 0]');
        obj.p = p_list(i);
    
        obj_list = [obj_list, obj];
    
        pc = get_pcost(obj.get_phi(), p_list(i), q, con_nor, con_dist);
        c_list = [c_list, pc];
    end
    
    [minc, minci] = min(c_list);
    obj = obj_list(minci);
end

p = obj.p;
ep = atanh((p - 2)*2/(obj.maxp-2) - 1);
obj.ep = p;
obj.dpdep = (obj.maxp-2)/2 * (1 - tanh(ep)^2);

end



%% Shaping cost
function [cost, grad] = shape_cost(var, q, con_nor_, con_dist_)
    tmpobj = DSF3d(var');

    Rot = quat2rotm(reshape(q(4:7), 1, [])); % Rotation
    trans = reshape(q(1:3), [3, 1]); % Translation

    if nargout > 1
        h_res = tmpobj.support(q, con_nor_, true, false, false, false);
        h = h_res.h;
        dhdphi = h_res.dhdphi;

        center_dist = con_nor_' * (Rot*tmpobj.center + trans); % (D,1)

        cost = h + center_dist - con_dist_';
        dcddphi = con_nor_' * Rot * tmpobj.dcdphi;
        grad = dhdphi + dcddphi;
    else
        h_res = tmpobj.support(q, con_nor_, false, false, false, false);

        center_dist = con_nor_' * (Rot*tmpobj.center + trans); % (D,1)

        cost = h_res.h + center_dist - con_dist_';
    end
end



function [cost, grad] = shape_cost_p(var, p, q, con_nor_, con_dist_)
    tmpobj = DSF3d([var; 0]');
    tmpobj.p = p;

    Rot = quat2rotm(reshape(q(4:7), 1, [])); % Rotation
    trans = reshape(q(1:3), [3, 1]); % Translation

    if nargout > 1
        h_res = tmpobj.support(q, con_nor_, true, false, false, false);
        h = h_res.h;
        dhdphi = h_res.dhdphi;

        center_dist = con_nor_' * (Rot*tmpobj.center + trans); % (D,1)

        cost = h + center_dist - con_dist_';
        dcddphi = con_nor_' * Rot * tmpobj.dcdphi;
        grad = dhdphi + dcddphi;
        grad = grad(:, 1:end-1);
    else
        h_res = tmpobj.support(q, con_nor_, false, false, false, false);

        center_dist = con_nor_' * (Rot*tmpobj.center + trans); % (D,1)

        cost = h_res.h + center_dist - con_dist_';
    end
end

function cost = get_pcost(var, p, q, con_nor_, con_dist_)
    tmpobj = DSF3d(var);
    tmpobj.p = p;

    Rot = quat2rotm(reshape(q(4:7), 1, [])); % Rotation
    trans = reshape(q(1:3), [3, 1]); % Translation

    h_res = tmpobj.support(q, con_nor_, false, false, false, false);

    center_dist = con_nor_' * (Rot*tmpobj.center + trans); % (D,1)

    cost = h_res.h + center_dist - con_dist_';

    cost = sum(cost.^2, "all");
end