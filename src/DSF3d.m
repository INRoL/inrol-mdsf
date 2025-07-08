classdef DSF3d
    properties
        v;
        vvt;
        loc_v; % (loc_v = v - center)

        alpha;
        center;
        dcdphi;

        ep; % (-inf~+inf)
        p; % (2~maxp)
        dpdep;

        maxp;
        Nv;
    end
    
    methods (Static)
    end

    methods
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Class Constructor
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Convex Constructor
        function obj = DSF3d(varargin)
            % varargin = {v, alpha, ep} or {phi}
            % v:(3,V)
            % alpha:(1,V)
            % center:(3,1)
            % dcdphi:(3,4*V+1)

            if nargin == 1
                phi = varargin{1};
                V = round((length(phi)-1)/4);
                v = reshape(phi(1:3*V), [3, V]);
                alpha = reshape(phi(3*V+1:4*V), [1, V]);
                ep = phi(end);
            elseif nargin == 3
                v = varargin{1};
                alpha = varargin{2};
                ep = varargin{3};
                V = size(v,2);
            else
                error('[Error] Wrong constructor input for DSF3d');
            end
            
            obj.Nv = V;
            obj.maxp = 40; % pre-determined

            obj.ep = ep;
            obj.p = (tanh(ep) + 1)*(obj.maxp-2)/2 + 2;
            obj.dpdep = (obj.maxp-2)/2 * (1 - tanh(ep)^2);
            
            obj.v = v; % (3,V)
            
            %----- interpolation
            obj.alpha = alpha;

            i_ub = 0.8; % upper bound
            a = (V-1)/(i_ub*V-1);
            k = exp(alpha);
            k = k / sum(k);
            interp = k/a + (a-1)/(a*V);
            dida = 1/a * (diag(k) - k'*k);

            %----- center
            obj.center = sum(v.*interp, 2); % (3,1)

            dcdv = reshape(eye(3), [3, 3, 1]) .* reshape(interp, [1, 1, V]);
            dcdv = reshape(dcdv, [3, 3*V]);

            dcda = v * dida;

            obj.dcdphi = [dcdv, dcda, zeros(3, 1)];
            
            %----- vvt in convex local frame
            obj.loc_v = v - obj.center;
            obj.vvt = reshape(obj.loc_v, [3, 1, V]) .* reshape(obj.loc_v, [1, 3, V]); % (3,3,V)
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Utilities
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Get shape vector phi
        function phi = get_phi(obj)
            phi = [reshape(obj.v, 1, []), obj.alpha, obj.ep];
        end
        %% Translation
        function obj = translate(obj, t)
            obj.v = obj.v + t;
            obj.center = obj.center + t;
        end
        %% Rotation
        function obj = rotate(obj, R)
            obj.v = R*obj.v;
            obj.center = R*obj.center;
            obj.loc_v = R*obj.loc_v;
        end
        %% Scaling
        function obj = scale(obj, scale)
            obj.v = obj.v * scale;
            obj.center = obj.center * scale;
            obj.dcdphi(:, (3*obj.Nv+1):(4*obj.Nv)) = scale * obj.dcdphi(:, (3*obj.Nv+1):(4*obj.Nv));
            obj.loc_v = obj.loc_v * scale;
            obj.vvt = obj.vvt * scale*scale;
        end
        %% Quaternion Transformation
        function obj = quat_transform(obj, q)
            Rot = quat2rotm(reshape(q(4:7), 1, [])); % Rotation
            trans = reshape(q(1:3), [3, 1]); % Translation
            obj.v = Rot*obj.v + trans;
            obj.center = Rot*obj.center + trans;
            obj.loc_v = Rot*obj.loc_v;
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Support Function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Support function
        function [h_res, s_res] = support(obj, q, x, get_dhdphi, get_s_rel, get_dsdx, get_dsdphi)
            % q = (1,7) : translation(3) + quaternion(4)
            % x = (3,D) : directions
            % get_dhdphi : true: return dhdphi
            % get_s_rel : true: return relative s (s-center) / false: s in global frame
            % get_dsdx : true: return dsdx
            % get_dsdphi : true: return dsdphi
            % h_res = struct[h, dhdphi] (D,1), (D,4*V+1)
            % s_res = struct[s, dsdx, dsdphi] (D,3), (D,3,3), (D,3,4*V+1)
            % dsdphi : always in convex local frame

            D = size(x, 2); % (3,D)
            V = obj.Nv;

            R01 = quat2rotm(reshape(q(4:7), 1, []));
            is_RI = false; % identity rotation
            if R01(1,1) == 1 && R01(2,2) == 1 && R01(3,3) == 1
                is_RI = true;
            end

            x1 = R01'*x;
            
            z = x1' * obj.loc_v; % (D,V)
            
            z = max(z, 0);
            maxz = max(z, [], 2); % (D,1)
            maxz = max(maxz, 1e-30);
            z = z ./ maxz;
            
            zp_2 = z .^ (obj.p - 2);
            zp_1 = zp_2 .* z; % (D,V)
            zp = zp_1 .* z;
            sum_zp = sum(zp, 2); % (D,1)
            
            h_2 = sum_zp .^ (1/obj.p - 2); % (D,1)
            h_1 = h_2 .* sum_zp;
            h = h_1 .* sum_zp;

            h_res.h = h .* maxz; % (D,1)
            
            if get_dhdphi == true
                dhdz = zp_1 .* h_1; % (D,V)
                dhdz = reshape(dhdz, [D, 1, V]);
                dhdlv = dhdz .* reshape(x1', [D, 3, 1]); % (D,3,V)
                dhdphi1 = [reshape(dhdlv, [D,3*V]), zeros(D,V)]; % (D,4*V) for v in v_loc

                sum_dhdlv = sum(dhdlv, 3); % (D,3,1)
                dhdphi2 = -reshape(sum_dhdlv, [D,3]) * obj.dcdphi(:, 1:4*V); % (D,4*V) for c in v_loc
                
                dhdphi = dhdphi1 + dhdphi2;
                
                z_forlog = z;
                z_forlog(z == 0) = 1;
                sum_zp_forlog = sum_zp;
                sum_zp_forlog(sum_zp == 0) = 1;
                log_z = log(z_forlog);
                log_sum_zp = log(sum_zp_forlog);
                
                sum_zplog = sum(zp .* log_z, 2); % (D,1)
                dhdp = h_1 .* (sum_zplog - log_sum_zp.*sum_zp/obj.p) / obj.p;
                dhdp = dhdp .* maxz; % (D,1)
                dhdep = dhdp * obj.dpdep;

                dhdphi = [dhdphi, dhdep]; % (D,4*V+1)

                h_res.dhdphi = dhdphi;
            end
 
            if nargout > 1
                zp_1_v = zp_1*obj.loc_v'; % (D,3)
                s = h_1 .* zp_1_v; % (D,3)
                if get_s_rel == false
                    s = s + obj.center'; % (D,3)
                    s0 = reshape(q(1:3), [1, 3]) + s*R01'; % global s ((R01*s')' = s*R01')
                else
                    s0 = s*R01'; % R01*(s - center)
                end
                s_res.s = s0;
        
                if get_dsdx == true
                    dsdx1 = reshape(zp_2, [1,D,1,V]) .* reshape(obj.vvt, [3,1,3,V]);
                    dsdx1 = reshape(sum_zp, [1,D,1]) .* reshape(sum(dsdx1, 4), [3,D,3]); % (3,D,3)

                    dsdx2 = reshape(zp_1_v', [3,D,1]) .* reshape(zp_1_v, [1,D,3]); % (3,D,3)
                    dsdx = dsdx1 - dsdx2;
                    dsdx = (obj.p-1) * reshape(h_2, [1,D,1]) .* dsdx; % (3,D,3)
                    dsdx = dsdx ./ reshape(maxz, [1,D,1]);
                    
                    if is_RI == false
                        dsdx0 = reshape(reshape(dsdx, [3*D, 3]) * R01', [3, D, 3]);
                        ds0dx0 = reshape(R01 * reshape(dsdx0, [3, 3*D]), [3, D, 3]);
                        ds0dx0 = permute(ds0dx0, [2, 1, 3]); % (D,3,3)
                        s_res.dsdx = ds0dx0;
                    else
                        dsdx = permute(dsdx, [2, 1, 3]); % (D,3,3)
                        s_res.dsdx = dsdx;
                    end
                end
        
                if get_dsdphi == true
                    %----- for local s
                    dsdlv = reshape(zp_2,[1,D,1,V]) .* reshape(obj.loc_v,[3,1,1,V]) .* reshape(sum_zp,[1,D,1,1]); % (3,D,1,V)
                    dsdlv = dsdlv - reshape(zp_1,[1,D,1,V]) .* reshape(zp_1_v',[3,D,1,1]); % (3,D,1,V)
                    dsdlv = (obj.p-1) * dsdlv .* reshape(x1', [1, D, 3, 1]); % (3,D,3,V)
                    dsdlv = dsdlv ./ reshape(maxz, [1, D, 1, 1]);
    
                    dsdlv2 = reshape(eye(3), [3, 1, 3, 1]) .* reshape(zp_1, [1, D, 1, V]) .* reshape(sum_zp,[1,D,1,1]);
                    dsdlv = dsdlv + dsdlv2;
                    dsdlv = reshape(h_2, [1, D, 1, 1]) .* dsdlv; % (3,D,3,V)
                    
                    dsdphi1 = cat(3, reshape(dsdlv, [3,D,3*V]), zeros(3,D,V)); % (3,D,4*V)

                    sum_dsdlv = sum(dsdlv, 4); % (3,D,3,1)
                    dsdphi2 = -reshape(sum_dsdlv, [3*D, 3])*obj.dcdphi(:, 1:4*V); % (3*D,4*V)
                    dsdphi2 = reshape(dsdphi2, [3, D, 4*V]); % (3,D,4*V)

                    dsdphi = dsdphi1 + dsdphi2;
                    
                    if get_dhdphi == false
                        z_forlog = z;
                        z_forlog(z == 0) = 1;
                        sum_zp_forlog = sum_zp;
                        sum_zp_forlog(sum_zp == 0) = 1;
                        log_z = log(z_forlog);
                        log_sum_zp = log(sum_zp_forlog);

                        sum_zplog = sum(zp .* log_z, 2); % (D,1)
                    end

                    dsdp = (1/obj.p-1)*sum_zplog;
                    dsdp = dsdp - log_sum_zp/(obj.p^2) .* sum_zp; % (D,1)
                    dsdp = dsdp .* zp_1_v; % (D,3)
    
                    dsdp = dsdp + ((zp_1.*log_z) * obj.loc_v') .* sum_zp;
    
                    dsdp = h_2 .* dsdp; % (D,3)
                    dsdep = dsdp * obj.dpdep;

                    dsdphi = cat(3, dsdphi, reshape(dsdep', [3,D,1])); % (3,D,4*V+1)
                    
                    if is_RI == false
                        ds0dphi = reshape(R01 * reshape(dsdphi, [3, D*(4*V+1)]), [3, D, 4*V+1]);
                        ds0dphi = permute(ds0dphi, [2, 1, 3]); % (D,3,4*V+1)
                        s_res.dsdphi = ds0dphi;
                    else
                        dsdphi = permute(dsdphi, [2, 1, 3]); % (D,3,4*V+1)
                        s_res.dsdphi = dsdphi;
                    end
                end
            end
        end

        %% Approximate distances
        function [dist, didx, dddphi] = approx_dist(obj, q, x, points)
            Rot = quat2rotm(reshape(q(4:7), 1, [])); % Rotation
            trans = reshape(q(1:3), [3, 1]); % Translation

            if nargout > 2
                h_res = obj.support(q, x, true, false, false, false);
                h = h_res.h;
                dhdphi = h_res.dhdphi;

                local_points = points - (Rot*obj.center + trans); % (3,P)

                p_dist = local_points' * x; % (P,D)
                dist = p_dist - h'; % (P,D)
                
                [dist, didx] = max(dist, [], 2); % (P,1)
                
                max_x = x(:, didx); % (3,P)
                max_dhdphi = dhdphi(didx, :); % (P,var)

                dddphi = -max_x'*(Rot*obj.dcdphi) - max_dhdphi;
            else
                h_res = obj.support(q, x, false, false, false, false);
                local_points = points - (Rot*obj.center + trans);
                p_dist = local_points' * x;
                dist = p_dist - h_res.h';
                [dist, didx] = max(dist, [], 2); % (P,1)
            end
        end

        %% Single direction
        function [h_res, s_res] = support_single(obj, q, x, get_dhdphi, get_s_rel, get_dsdx, get_dsdphi)
            % More efficient way
            R01 = quat2rotm(q(4:7));
            x1 = R01'*x;
            z = x1' * obj.loc_v; % (1,V)
            z = max(z, 0);
            max_z = max(z);
            z = z/max_z;
            
            zp_2 = z.^(obj.p-2);
            zp_1 = zp_2.*z;
            zp = zp_1.*z;

            sum_zp = sum(zp);
            h_2 = sum_zp^(1/obj.p - 2);
            h_1 = h_2*sum_zp;
            h = h_1*sum_zp;
            
            %----- h
            h_res.h = h * max_z;

            if get_dhdphi == true
                dhdv = h_1*zp_1 .* x1; % (3,V)
                dhda = -sum(dhdv, 2)' * obj.dcdphi; % (1,4*V+1)
                
                z_forlog = z;
                z_forlog(z == 0) = 1;
                sum_zp_forlog = sum_zp;
                sum_zp_forlog(sum_zp == 0) = 1;
                log_z = log(z_forlog);
                log_sum_zp = log(sum_zp_forlog);
                sum_zplog = sum(zp.*log_z);

                dhdp = h_1*sum_zplog/obj.p - h*log_sum_zp/(obj.p^2);
                dhdep = dhdp * obj.dpdep * max_z;

                dhdphi = [reshape(dhdv, 1, []), zeros(1, obj.Nv+1)];
                dhdphi = dhdphi + dhda;
                dhdphi(end) = dhdep;

                h_res.dhdphi = dhdphi;
            end

            if nargout > 1
                sum_zp_1_v = obj.loc_v * zp_1'; % (3,1)
                s = h_1 * sum_zp_1_v; % (3,1)
                if get_s_rel == false
                    s = s + obj.center;
                    s = q(1:3)' + R01*s;
                else
                    s = R01*s;
                end
                s_res.s = s;

                if get_dsdx == true
                    dsdx = h_1 * sum(reshape(zp_2, [1,1,obj.Nv]) .* obj.vvt, 3);
                    dsdx = reshape(dsdx, [3,3]) - h_2*(sum_zp_1_v * sum_zp_1_v');
                    dsdx = (obj.p-1) * R01 * dsdx * R01';
                    s_res.dsdx = dsdx / max_z;
                end

                if get_dsdphi == true
                    dsdlv = h_1*zp_2.*obj.loc_v; % (3,V)
                    dsdlv = dsdlv - h_2*zp_1.*sum_zp_1_v; % (3,V)
                    dsdlv = (obj.p-1)*dsdlv;
                    dsdlv = reshape(dsdlv, [3, 1, obj.Nv]) .* reshape(x1, [1,3,1]); % (3,3,V)
                    dsdlv = dsdlv / max_z + h_1*reshape(eye(3),[3,3,1]).*reshape(zp_1, [1,1,obj.Nv]);
                    dsdv = [reshape(dsdlv, [3, 3*obj.Nv]), zeros(3, obj.Nv+1)]; % (3,4*V+1)
                    
                    dsda = -reshape(sum(dsdlv, 3), [3,3]) * obj.dcdphi; % (3,4*V+1)
                    
                    if get_dhdphi == false
                        z_forlog = z;
                        z_forlog(z == 0) = 1;
                        sum_zp_forlog = sum_zp;
                        sum_zp_forlog(sum_zp == 0) = 1;
                        log_z = log(z_forlog);
                        log_sum_zp = log(sum_zp_forlog);
                        sum_zplog = sum(zp.*log_z);
                    end
                    dsdp = ((1/obj.p-1)*h_2*sum_zplog - h_1*log_sum_zp/(obj.p^2)) * sum_zp_1_v; % (3,1)
                    dsdp = dsdp + h_1* obj.loc_v*(zp_1.*log_z)';
                    dsdep = dsdp * obj.dpdep;

                    dsdphi = dsdv + dsda;
                    dsdphi(:, end) = dsdep;
                    dsdphi = R01*dsdphi;
                    s_res.dsdphi = dsdphi;
                end
            end
        end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Visualizations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Visualization
        function draw(obj, q, surfdir, n1, n2, color, alpha, edgealpha, show_v, show_t)
            Rot = quat2rotm(reshape(q(4:7), 1, [])); % Rotation
            trans = reshape(q(1:3), [3, 1]); % Translation
            %----- show vertices
            if show_v == true
                tmpv = Rot*obj.v + trans;
                scatter3(tmpv(1, :), tmpv(2, :), tmpv(3, :), 300, "r", "filled");
            end
            
            %----- show center
            if show_t == true
                c = Rot*obj.center + trans;
                scatter3(c(1), c(2), c(3), "g", "filled");
            end
            
            %----- show surface
            [unused, s_res] = obj.support(q, surfdir, false, false, false, false);

            s = reshape(s_res.s, [n1, n2, 3]);
            surf(s(:, :, 1), s(:, :, 2), s(:, :, 3),'Edgecolor', "none", 'Facealpha',alpha,'Edgealpha',edgealpha,'facecolor',color);
        end
    end
end