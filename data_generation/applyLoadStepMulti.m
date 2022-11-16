function [flag_solv,DISP,sigma_bar,bar_to_remove,K,C,res_disp_fin, ind] = applyLoadStepMulti(A_bar,...
    E_conc,Load_steps,force_vec,K,partion_1,tol_iter,max_iter,kk,...
    Nodes,Con_bar,flag_solv,DISP,thetaX,thetaY,thetaZ,partion_2,...
    tensile_strength,bar_to_remove,Ig_bar,Jg_bar,L_bar,K_tri,K_quad,radius)

E_conc_vec = ones(size(A_bar))*E_conc;
E_diag_con = E_conc_vec;
total_load = 0;
sigma_bar = 0;
ii = 1;
p = 1;
C = {};
res_disp_fin = 0;
ind = 0;
tic
while p/Load_steps < 1
%     disp("load step " + num2str(ii))
    n = 1;
    force_vec_new = (force_vec./Load_steps) + total_load;
    no_bar_to_remove = 1;
    K_old = K;
%     disp("percent complete: " + num2str(100 * p/Load_steps) + "%")
    while no_bar_to_remove > 0
        % MZ check to make sure # of load steps is reasonable
        if length(bar_to_remove{ii,kk}) > 100
            Load_steps = Load_steps + 10;
            bar_to_remove{ii,kk} = [];
            flag_solv{ii,kk} = [];
            if kk == 1 || length(bar_to_remove) < Load_steps
                old_length = length(bar_to_remove);
                bar_to_remove_new = cell(Load_steps,3);
                bar_to_remove_new(1:old_length,1:kk) = bar_to_remove(:,1:kk);
                bar_to_remove = bar_to_remove_new;
                flag_solv_new = cell(Load_steps,3);
                flag_solv_new(1:old_length,1:kk) = flag_solv(:,1:kk);
                flag_solv = flag_solv_new;
            end
            K = K_old;
            p = (p+10) - 1;
            break
        end
        %wu ignore the support dof
%         disp("iteration " + num2str(n))
        K_solve = K;
        K_solve(partion_1,:) = [];
        K_solve(:, partion_1) = [];
        
        force_vec_solv = force_vec_new;
        force_vec_solv(partion_1,:) = [];
        
        %wu solve kx=f
        % %                 [res_disp, flag_solv{ii,kk}(n,1),relres,iter,resvec] = cgs(K_solve,force_vec_solv,tol_iter,max_iter);
        %                 tic
        %                 K_solve=round(K_solve,5);
        %                 K_solve_MF = Multifrontal(K_solve);
        %                 res_disp2=K_solve_MF\force_vec_solv;
        %                 toc
        
        tic
        K_solve_dist = gpuArray(K_solve);
        force_vec_solv_dist = gpuArray(force_vec_solv);
        [res_disp, flag_solv{ii,kk}(n,1)] = cgs(K_solve_dist,force_vec_solv_dist,tol_iter,max_iter);
%         res_disp = gather(res_disp);
%         flag_solv{ii,kk}(n,1) = gather(flag_solv{ii,kk}(n,1));
%         delete(gcp('nocreate'));
        toc;
        
        
        %                 flag_solv{ii,kk}(n,1)
        %                 relres
        %                 numel_resvec=numel(resvec)
        %                 my_relres=norm(K_solve*res_disp-force_vec_solv)/norm(force_vec_solv)
        
        %wu displacement of nodes
        res_disp_fin = zeros(3*size(Nodes,1),1);
        res_disp_fin(partion_2,:) = res_disp;
        
        %wu ??why magnify??
        DISP{ii,kk} = reshape(res_disp_fin,3,[])';
        
        %wu displacement of a bar element
        u_bar(1,:) = res_disp_fin(Con_bar(:,1)*3-2,1);
        u_bar(2,:) = res_disp_fin(Con_bar(:,1)*3-1,1);
        u_bar(3,:) = res_disp_fin(Con_bar(:,1)*3,1);
        u_bar(4,:) = res_disp_fin(Con_bar(:,2)*3-2,1);
        u_bar(5,:) = res_disp_fin(Con_bar(:,2)*3-1,1);
        u_bar(6,:) = res_disp_fin(Con_bar(:,2)*3,1);
        
        %MZ find stress of bar elements
        sigma_bar = zeros(size(Con_bar,1),1);
        C = cell(size(Con_bar,1),1);
        for i = 1:size(Con_bar,1)
            C_temp = [-thetaX(i,1) -thetaY(i,1) -thetaZ(i,1) thetaX(i,1) thetaY(i,1) thetaZ(i,1)];
            sigma_bar(i,1) = (E_conc_vec(i,1)/L_bar(i,1))*C_temp*u_bar(:,i);
            C(i,1) = {C_temp};
        end
        for j = 1:1
            disp("Element Removed: " + num2str(j))
            max_sigma_val = max(sigma_bar);
            max_sigma_ind = find(sigma_bar == max_sigma_val);
            sigma_bar(max_sigma_ind,1) = 0;
        
        %wu find max stress and remove
            if max_sigma_val > tensile_strength(max_sigma_ind,1)
                bar_to_remove{ii,kk}(1,n) = max_sigma_ind(1,1);
                E_conc_vec(bar_to_remove{ii,kk}(1,n),1) = 0;
                E_diag_con(bar_to_remove{ii,kk}(1,n),1) = 1e-6;
            
                const = E_conc_vec.*A_bar./L_bar;
                const_diag = E_diag_con.*A_bar./L_bar;
                kc_prime = [1 -1; -1 1];
                Kg_bar_m = cell(size(Con_bar,1),1); 
                for i = 1:size(Con_bar,1)
                    const_m = [const_diag(i,1) const(i,1) const(i,1) const(i,1) const(i,1) const(i,1);
                    const(i,1) const_diag(i,1) const(i,1) const(i,1) const(i,1) const(i,1);
                    const(i,1) const(i,1) const_diag(i,1) const(i,1) const(i,1) const(i,1);
                    const(i,1) const(i,1) const(i,1) const_diag(i,1) const(i,1) const(i,1);
                    const(i,1) const(i,1) const(i,1) const(i,1) const_diag(i,1) const(i,1);
                    const(i,1) const(i,1) const(i,1) const(i,1) const(i,1) const_diag(i,1)];
                    Tstar = [thetaX(i,1) 0; thetaY(i,1) 0; thetaZ(i,1) 0; 0 thetaX(i,1); 0 thetaY(i,1); 0 thetaZ(i,1)];
                    Kg_bar_m(i,1) = {const_m.*(Tstar*kc_prime*Tstar')};
                    Kg_bar_m{i,1} = reshape(Kg_bar_m{i,1}, [], 1);
                end
                Kg_bar_new = cat(1, Kg_bar_m{:});
            
                nq = size(Nodes,1);
                Kg_bar_new_sparse = sparse(Ig_bar(:),Jg_bar(:),Kg_bar_new(:),3*nq,3*nq);
                K = Kg_bar_new_sparse+K_tri+K_quad;
            
                no_bar_to_remove = 1;
                n = n+1;
                
                x = [Nodes(Con_bar(bar_to_remove{ii,kk}(1,n),1),1); Nodes(Con_bar(bar_to_remove{ii,kk}(1,n),2),1)];
                y = [Nodes(Con_bar(bar_to_remove{ii,kk}(1,n),1),2); Nodes(Con_bar(bar_to_remove{ii,kk}(1,n),2),2)];
                z = [Nodes(Con_bar(bar_to_remove{ii,kk}(1,n),1),3); Nodes(Con_bar(bar_to_remove{ii,kk}(1,n),2),3)];
                
                x_cent = x - radius;
                y_cent = y - radius;
                
                [~, rho, ~] = cart2pol(x_cent,y_cent,z);
                if max(rho) >= radius
                    no_bar_to_remove = 0;
                    ind = 1;
                    p = Load_steps + 1;
                    break
                end
            else
                no_bar_to_remove = 0;
                break
            end
        end
    end
    if ~(isempty(bar_to_remove{ii,kk})) || n==1
        total_load = force_vec_new;
        ii = ii + 1;
        p = 1 + p;
    end
end
