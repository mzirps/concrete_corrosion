% function crack_main(geometry_info)
clc
clear all
close all
%reference
%[1] MULTI-PHYSICAL AND MULTI-SCALE DETERIORATION MODELLING OF RE-INFORCED CONCRETE: MODELLING CORROSION-INDUCED CONCRETE DAMAGE
cd '/home/ubuntu/sherlock_data/GIBBON'
installGibbon
cd '/home/ubuntu/sherlock_data/dir'

%%
 %m
elemlength=.5; %m


pot_rebarR = [10, 12, 13, 14, 16]/2/10^3;

%MZ input constants for mesh info
% number of elements around circumference of rebar
No_ele_rebar_surf = 30;
No_ele_rebar_vert = 80;
% number of random nodes in concrete domain
conc_no_spacing_bound = 0.015;
no = 40;
No_nodes_x = no;
No_nodes_y = no;
No_nodes_z = no;
% number of random nodes in rebar domain
No_nodes_x_rebar = 10;
No_nodes_y_rebar = 10;
No_nodes_z_rebar = 10;

warning('off','all')
% asd
% all dimensions in [mm]
%MZ input origin
X_origin = 0;
Y_origin = 0; 
Z_origin = 0;

%MZ input material properties
% all dimensions in [MPa]
E_conc = 24300;
E_steel = 210000;
NU_steel = 0.2;
E_corr = 2000;
NU_corr = 0.2;

Load_steps = 30;
% fitting parameter for area calculation of bar elements  
fit_1 = 50;

% size of sub-cells for creation of random points, between 0.5 (ordered)
% and 1 (random)
sub_cell_size = 0.6;
Magnification_deformation = 10;

% input parameters penetration of corrosion products
CAR_0 = 0.14;
CAR_max = 0.28;
w_c = 0.74;
n_non = 1.3;
delta_T = 0.7;

timeSteps=1;


for i = 7:9450
    t = 1;
    if i == 7
        t = 8;
    end
    while t < 10
        rebarR = pot_rebarR(1,randi(5));
        cover = 0.06 + rand*0.04;
        radius = cover + rebarR;
        tensile_strength = 2 + rand*3;
        w_c = .4 + rand*.2;
        
        longRebar = cell(2,1);
        longRebar{1} = [radius radius 0 radius radius elemlength];
        longRebar{2} = [1 2*rebarR 0 0 radius radius radius radius];
        
        num_longRebar=sum(longRebar{2,1}(:,1));

        [deltax,deltay,deltaz,rebarR,rebarCenter] = getRebarDim(longRebar,num_longRebar);
        
        l_A = deltaz;
        
        if w_c < 0.42
            m = w_c/0.42;
        else
            m = 0.65;
        end
        coeff_por = (w_c-0.36*m)/(w_c+0.32);
        V_cp_min = coeff_por*deltaz.*pi()*1000.*((rebarR*1000).^2 - ((rebarR*1000) - 0.11).^2); %mm3
        V_cp_max = coeff_por*deltaz.*pi()*1000.*((rebarR*1000).^2 - ((rebarR*1000) - 0.45).^2); %mm3

        filename = "/home/ubuntu/corrosion/Corrosion_simulation_" + num2str(i) + "_timeStep_" + num2str(t) + ".txt";
        exp_data = importdata(filename);
        exp_data = cell2num(exp_data(4:end,:));
        
        % MZ new experimental data
        exp_corr_layer=cell(num_longRebar,timeSteps);
        exp_theta_corr_layer=cell(num_longRebar,1);
        exp_z_corr_layer = cell(num_longRebar,1);

        exp_corr_layer{1,1} = [exp_data(:,2); exp_data(:,2)];
        if t ~= 1
            %mult by 5/10^6 for trials 1-6 1/10^5 for 7-
            exp_corr_layer{1,1} = [exp_data(:,2); exp_data(:,2)]*(1/10^5);
        end
        exp_theta_corr_layer{1,1} = [0; 359];
        exp_z_corr_layer{1,1} = exp_data(:,1);
        
        % create random points in sub-cells
        bar_to_remove = cell(Load_steps,timeSteps);
        DISP = cell(Load_steps,timeSteps);
        flag_solv = cell(Load_steps,timeSteps);
        
        non_uni_corr_layer=cell(num_longRebar,timeSteps);
        
        max_iter = 5e4;
        tol_iter = 1e-3;
        
        [bar_to_remove, theta, z, ind] = Crack3D(num_longRebar,rebarR,No_ele_rebar_surf,No_ele_rebar_vert,...
            deltax,deltay,deltaz,rebarCenter,exp_z_corr_layer,exp_theta_corr_layer,...
            delta_T,l_A,V_cp_min,V_cp_max,CAR_0,CAR_max,coeff_por,n_non,...
            exp_corr_layer,non_uni_corr_layer,radius,elemlength,conc_no_spacing_bound,...
            X_origin,Y_origin,Z_origin,No_nodes_x,No_nodes_y,No_nodes_z,sub_cell_size,...
            No_nodes_x_rebar,No_nodes_y_rebar,No_nodes_z_rebar,fit_1,E_conc,tensile_strength,...
            E_steel,NU_steel,E_corr,NU_corr,Load_steps,tol_iter,max_iter,flag_solv,DISP,...
            bar_to_remove);
        name = "Data_outputs/output_" + num2str(i) + "_" + num2str(t);
        parsave(name,rebarR,cover,tensile_strength,w_c,bar_to_remove,theta,z,ind)
        t = t+1;
        
    end
end

