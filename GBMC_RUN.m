%  ----------------------------------------------------------------
% RESEARCH ONLY LICENSE
% Copyright © 2018 North Carolina State University.
% All rights reserved.
% Redistribution and use in source and binary forms, with or without 
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions and use are permitted for internal research purposes only, 
% and commercial use is strictly prohibited under this license. Inquiries 
% regarding commercial use should be directed to the Office of Technology 
% Transfer at North Carolina State University, 919‐215‐7199, https://research.
% ncsu.edu/ott/contact‐us/techtransfer@ncsu.edu.
% 
% 2. Commercial use means the sale, lease, export, transfer, conveyance or other 
% distribution to a third party for financial gain, income generation or other 
% commercial purposes of any kind, whether direct or indirect. Commercial use 
% also means providing a service to a third party for financial gain, income 
% generation or other commercial purposes of any kind, whether direct or 
% indirect.
% 
% 3. Redistributions of source code must retain the above copyright notice, this 
% list of conditions and the following disclaimer.
% 
% 4. Redistributions in binary form must reproduce the above copyright notice, 
% this list of conditions and the following disclaimer in the documentation and/
% or other materials provided with the distribution.
% 
% 5. The names “North Carolina State University”, “NCSU” and any trade‐name,
% personal name, trademark, trade device, service mark, symbol, image, icon, or 
% any abbreviation, contraction or simulation thereof owned by North Carolina 
% State University must not be used to endorse or promote products derived from 
% this software without prior written permission. For written permission, please
% contact trademarks@ncsu.edu.
% 
% Disclaimer: THIS SOFTWARE IS PROVIDED “AS IS” AND ANY EXPRESSED OR IMPLIED 
% WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
% MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO 
% EVENT SHALL THE APACHE SOFTWARE FOUNDATION OR ITS CONTRIBUTORS BE LIABLE FOR 
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES 
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; 
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON 
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (
% INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% -------------------------------------------------------------------------
% Authors: Arash Dehghan Banadaki (adehgha@ncsu.edu)
%          Srikanth Patala        (spatala@ncsu.edu)
% -------------------------------------------------------------------------
% This software is the implementation of the Monte-Carlo (random insertion-
% removal) technique for Grain Boundary simulation that has been presented in the 
% following article. Please consider citing this article if you have used this 
% software in your work. 
% [Banadaki A.D., Tschopp M.A., Patala S., An efficient Monte Carlo algorithm 
% for determining the minimum energy structures of metallic grain boundaries,
% Computational Materials Science, v.155, p.466-475]
% https://doi.org/10.1016/j.commatsci.2018.09.017
% -------------------------------------------------------------------------

function GBMC_RUN(gb_main_folder)

% ===============S I M U L A T I O N======================
% +===============C O N S T A N T S======================

lammps_exe_path = '/usr/bin/lmp_daily';
max_steps = 5000;
insertion_freq = 2;
global lammps_exe_path;
T = 500; % K
% --------for ni simulations ----------
% elem = 'ni';
% csm_crit = 0.1;
% lat_type = 'fcc';
% lat_param = 3.52;
% --------for al simulations ----------
elem = 'al';
csm_crit = 0.1;
lat_type = 'fcc';
lat_param = 4.05;


% ========================================================
K_boltz = 1.3806485279e-23; % J.K^-1
event_num = 1;

if isdir(gb_main_folder)
        gb_structure_file = dir([gb_main_folder '/dump.*']);
        dump_file        = fullfile(gb_main_folder, ['dump_1']);
        copyfile(fullfile(gb_main_folder, gb_structure_file.name), dump_file)
        
        eng_current = 1e5; %some large number
        eng_new = 0;
        cnt = 0;
        new_case_flag = true;
        all_acc_move_engs = zeros(max_steps, 2);
        extra_info = zeros(max_steps, 11);
        eng_file_addr = fullfile(gb_main_folder, 'all_acc_move_engs.mat');
        extra_file_addr = fullfile(gb_main_folder, 'extra_info.mat');
        header_info = ['eng_new ','new_case_flag ','mode_num ',...
                'perfect_cryst_depth_l ', 'min_perfect_cryst_depth_l ',...
                'perfect_cryst_depth_u ', 'min_perfect_cryst_depth_u ',...
                'mid_y_u ', 'max_y_csm ', 'min_y_csm ', 'mid_y_l '];
        gb_area_uc = reag_gb_area(fullfile(gb_main_folder, 'gb_area'));
        for cnt1 = 1:max_steps
            move_rand_num = rand();
            
            if cnt1 == 1
                mode_num = 0;
            else
                if move_rand_num < 0.50
                    mode_num = 6;
                else
                    mode_num = 1;
                end
            end
            
            
            tmp_file_name = fullfile(gb_main_folder, ['dump_' num2str(cnt1)]);
            %     if ~exist(tmp_file_name, 'file') == 2
            if new_case_flag
                dump_file = tmp_file_name;
                [MC_data] = MC1(dump_file, elem, lat_param, csm_crit, lat_type);
                if cnt1 == 1
                    [gb_attr] = find_gb_attributes(MC_data.Coords, csm_crit);
                    [l_p_go] = l_p_go_finder(MC_data, gb_attr, csm_crit, lat_param, lat_type);
                end
            end
            
            if mode_num ~= 7
                atom_remove_insertion(MC_data, gb_main_folder, mode_num, cnt1);
            else
                in_plane_shift(MC_data, gb_main_folder, csm_crit, cnt1);
                tic;
                move_back(lat_param, gb_main_folder, cnt1);
                toc;
            end
            
            if cnt1 ~= 1
                eng_current = eng_new;
            end
            gb_eng_file = fullfile(gb_main_folder, 'GB_Eng.dat');
            if exist(gb_eng_file, 'file') == 2
                [eng_new, ~, success_flag, extra_data] = eng_read(gb_eng_file);
                delta_eng = eng_new - eng_current;
                if delta_eng < 0 && success_flag
                    cnt = cnt+1;
                    new_case_flag = true;
                    all_acc_move_engs(cnt, :) = [eng_new, cnt1];
                elseif delta_eng >= 0 && success_flag
                    E_acc_prb = exp(-(delta_eng * gb_area_uc * 1e-20)/(K_boltz * 1e3 * T));
                    rand_k = rand();
                    if rand_k < E_acc_prb
                        cnt = cnt + 1;
                        new_case_flag = true;
                        all_acc_move_engs(cnt, :) = [eng_new, cnt1];
                    else
                        new_case_flag = false;
                        delete(fullfile(gb_main_folder, ['dump_' num2str(cnt1+1)]));
                    end
                else
                    %                     break;
                    add_single_crystal_block(elem, lat_param, csm_crit * 0.01, gb_main_folder, l_p_go, cnt1);
                    if exist(gb_eng_file, 'file') == 2
                        [eng_new, ~, success_flag, extra_data] = eng_read(gb_eng_file);
                        delta_eng = eng_new - eng_current;
                        if delta_eng < 0 && success_flag
                            cnt = cnt+1;
                            new_case_flag = true;
                            all_acc_move_engs(cnt, :) = [eng_new, cnt1];
                        elseif delta_eng >= 0 && success_flag
                            E_acc_prb = exp(-(delta_eng * gb_area_uc * 1e-20)/(K_boltz * 1e3 * T));
                            rand_k = rand();
                            if rand_k < E_acc_prb
                                cnt = cnt + 1;
                                new_case_flag = true;
                                all_acc_move_engs(cnt, :) = [eng_new, cnt1];
                            else
                                new_case_flag = true;
                            end
                        end
                    end
                end
            else
                break;
            end
            extra_info(cnt1,:) = [eng_new, new_case_flag, mode_num, extra_data];
            save(extra_file_addr,'extra_info','header_info');
            save(eng_file_addr,'all_acc_move_engs');
        end
        intrmd_steps_folder = fullfile(gb_main_folder,'intermediate_steps');
        mkdir(intrmd_steps_folder);
        movefile(fullfile(gb_main_folder,'dump_*.t'), intrmd_steps_folder);
        if cnt < max_steps && mod(cnt, 10)==0
            all_acc_move_engs(cnt+1:end,:) = [];
            save(eng_file_addr,'all_acc_move_engs');
            
            extra_info(cnt+1:end,:) = [];
            save(extra_file_addr,'extra_info','header_info');
        end
    end
    accepted_dumps = fullfile(gb_main_folder,'accepted_steps');
    mkdir(accepted_dumps);
    movefile(fullfile(gb_main_folder,'dump_*'), accepted_dumps);
    display(gb_main_folder);
    cd(gb_main_folder);
    
    cmnd_zip_1 = ['tar -zcf ', ['intermediate_steps.tar.gz '] './intermediate_steps'];
    status = system(cmnd_zip_1);
    
    cmnd_zip_2 = ['tar -zcf ', ['accepted_steps.tar.gz '] './accepted_steps'];
    status = system(cmnd_zip_2);
    if exist('./intermediate_steps.tar.gz','file') == 2
        cmnd_rm_1 = ['rm -rf ', intrmd_steps_folder];
        status = system(cmnd_rm_1);
    end
    if exist('./accepted_steps.tar.gz','file') == 2
        cmnd_rm_2 = ['rm -rf ', accepted_dumps];
        status = system(cmnd_rm_2);
    end
%end
end
%--------------------------------------------------------------------------


function [eng_new, gb_area, success_flag, extra_data] = eng_read(gb_eng_file)
% gb_eng_file = fullfile(gb_main_folder, 'GB_Eng.dat');
fid = fopen(gb_eng_file, 'r');
line = fgetl(fid); fclose(fid);
data = str2num(line);
eng_new = data(1);
gb_area = data(2);
if data(3) ==1
    success_flag = true;
else
    success_flag = false;
end
extra_data = data(4:end);
delete(gb_eng_file);
end
%--------------------------------------------------------------------------


function [MC_data] = MC1(gb_dump_address, elem, lat_param, CSMcrt, lat_type)
% CSMcrt = 0.1;
if strcmpi(elem,'al')
    min_energy = -3.36000010077;
    lat_param = 4.05;
    lat_type = 'fcc';
elseif strcmpi(elem,'cu')
    min_energy = -3.53999996838;
    lat_param = 3.615;
    lat_type = 'fcc';
elseif strcmpi(elem,'ni')
    min_energy = -4.45000000527;
    lat_param = 3.52;
    lat_type = 'fcc';
elseif strcmpi(elem,'fe')
    min_energy = -4.1224310073005;
    lat_param = 2.855324;    
    lat_type = 'bcc';
elseif strcmpi(elem,'mo')
    min_energy = -6.8200023329166;
    lat_param = 3.1472;    
    lat_type = 'bcc';
elseif strcmpi(elem,'si')
    min_energy = -4.63041206421337;
    lat_param = 5.431;    
    lat_type = 'diamond';
end

if exist(gb_dump_address)
    [Coords, dump_xyz_lims] = f_read(gb_dump_address);
    [high_csm_gb_atoms_indx, num_gb_atoms] = find_number_of_gb_atoms(Coords, dump_xyz_lims, CSMcrt);
    [rm_prb] = find_removal_prb(Coords, high_csm_gb_atoms_indx, min_energy);
    [insertion_centers, insertion_prb] = find_insertion_prb(Coords, dump_xyz_lims, CSMcrt, lat_param, lat_type);
    MC_data.elem              = elem;
    MC_data.Coords            = Coords;
    MC_data.dump_xyz_lims     = dump_xyz_lims;
    MC_data.rm_prb            = rm_prb;
    MC_data.insertion_centers = insertion_centers;
    MC_data.insertion_prb     = insertion_prb;
end
end
%--------------------------------------------------------------------------


function [] = atom_remove_insertion(MC_data, gb_main_folder, mode_num, cnt1)

max_mode_insert = 1;

elem = MC_data.elem;
if strcmpi(elem,'al')
    min_energy = -3.36000010077;
    lat_param = 4.05;
    lat_type = 'fcc';
elseif strcmpi(elem,'cu')
    min_energy = -3.53999996838;
    lat_param = 3.615;
    lat_type = 'fcc';
elseif strcmpi(elem,'ni')
    min_energy = -4.45000000527;
    lat_param = 3.52;
    lat_type = 'fcc';
elseif strcmpi(elem,'fe')
    min_energy = -4.1224310073005;
    lat_param = 2.855324;    
    lat_type = 'bcc';
elseif strcmpi(elem,'mo')
    min_energy = -6.8200023329166;
    lat_param = 3.1472;    
    lat_type = 'bcc';
end
% ---------------

Coords = MC_data.Coords;
dump_xyz_lims = MC_data.dump_xyz_lims;

% ---------------
switch mode_num
    case 0
    case {1,2,3,4}
        for mod_indx = 1:mode_num
            [Coords] = remove_atom_w_prb(MC_data);
        end
    case 5
        max_mode_remove = 1;
        if mod(cnt1,2)~=0
            [Coords, success_flag] = insert_atom_w_prb(MC_data);
        else
            for mod_indx = 1:max_mode_remove
                [Coords] = remove_atom_w_prb(MC_data);
            end
        end
    case 6
        [Coords] = insert_atom_w_prb(MC_data);
end

dump_file = fullfile(gb_main_folder, ['dump_' num2str(cnt1)]);
lammps_data_file_name = [dump_file, '.t'];
gb_eng_path = fullfile(gb_main_folder,'GB_Eng.dat');
lammps_data_file_write(Coords(:, 1:5), dump_xyz_lims, lammps_data_file_name);
% ----------------------------------
length_cnt1 = length(num2str(cnt1));
next_dump = [dump_file(1:end-length_cnt1) num2str(cnt1+1)];
pot_file = ['pot_', elem];
%cmnd1= 'mpirun -np 10 lammps-daily -log none -var';
global lammps_exe_path;
cmnd1= ['mpirun ', lammps_exe_path, ' -log none -var'];
cmnd2= [' potname1 ', pot_file,' -var potname2 ', pot_file,' -var '];
%         cmnd3 = ['fname1 ',lammps_data_file_name,' -var fname2 ',copper_dump, ' < in.GB_convert_AtoB.txt'];
cmnd3 = ['fname1 ',lammps_data_file_name,' -var fname2 ',next_dump,...
         ' -var fname3 ',gb_eng_path, ' -var min_eng ', num2str(min_energy,15),...
         ' -var lat_type ', lat_type, ' < in.minimize'];
cmnd = [cmnd1 cmnd2 cmnd3];
system(cmnd);
end

%--------------------------------------------------------------------------
function [Coords] = remove_atom_w_prb(MC_data)
%     [max_prb, indx_prb] = max(MC_data.rm_prb);
non_zero_prb_indx = find(MC_data.rm_prb);
non_zero_prb = MC_data.rm_prb(non_zero_prb_indx);
selected_indx = choose_1_among_n_moves_w_weighted_prb(non_zero_prb);
atom_to_remove = non_zero_prb_indx(selected_indx);
Coords = MC_data.Coords;
Coords(atom_to_remove,:) = [];
end
%--------------------------------------------------------------------------


function [Coords_w_extra_atom] = insert_atom_w_prb(MC_data)

non_zero_prb = MC_data.insertion_prb;
selected_indx = choose_1_among_n_moves_w_weighted_prb(non_zero_prb);

Coords = MC_data.Coords;
sph_cent = MC_data.insertion_centers(selected_indx, :);
Coords_w_extra_atom = [Coords ;Coords(end,:)];
Coords_w_extra_atom(end,3:5) = sph_cent;
Coords_w_extra_atom(end,1) = max(Coords_w_extra_atom(:,1))+1;
end
%--------------------------------------------------------------------------


function [selected_indx] = choose_1_among_n_moves_w_weighted_prb(prb)

[sorted_prb, i_prb] = sort(prb);
cum_sum_prb         = cumsum(sorted_prb);

rand_number = rand();
selected_indx = i_prb(min(find(cum_sum_prb > rand_number)));

end
%--------------------------------------------------------------------------


function [Coords, dump_xyz_lims]=f_read(NameOfFile)
% ---------------
fid=fopen(NameOfFile,'r');
tmp=fgetl(fid);tmp=fgetl(fid);tmp=fgetl(fid);tmp=fgetl(fid);
NoOfLines=str2num(tmp);
tmp=fgetl(fid);
if length(tmp)>25
    x_lims=fscanf(fid,'%g %g %g',[1 3]);
    y_lims=fscanf(fid,'%g %g %g',[1 3]);
    z_lims=fscanf(fid,'%g %g %g',[1 3]);
    a_tri = abs(x_lims(2)- x_lims(1))- abs(y_lims(3));
    c_tri = abs(z_lims(2)- z_lims(1));
    GB_area = a_tri * c_tri;
    tmp=fgetl(fid);
else
    x_lims=fscanf(fid,'%g %g %g',[1 2]);
    y_lims=fscanf(fid,'%g %g %g',[1 2]);
    z_lims=fscanf(fid,'%g %g %g',[1 2]);
    GB_area = abs(x_lims(2)-x_lims(1)) * abs(z_lims(2)-z_lims(1));
end
dump_xyz_lims = [x_lims; y_lims; z_lims];
tmp=fgetl(fid);
Coords=fscanf(fid,'%g %g %g %g %g %g %g',[NoOfLines 7]);
Coords=(reshape(Coords,7,NoOfLines))';
fclose(fid);
end
%--------------------------------------------------------------------------


function [rm_prb] = find_removal_prb(Coords, high_csm_gb_atoms, min_energy)
rm_prb = zeros(size(Coords, 1), 1);

positive_eng_atoms = find(Coords(:, 7) >= min_energy);
pos_eng_high_csm_gb_atoms = intersect(positive_eng_atoms,high_csm_gb_atoms);
rm_prb(pos_eng_high_csm_gb_atoms) = (Coords(pos_eng_high_csm_gb_atoms, 7)-min_energy)/sum(Coords(pos_eng_high_csm_gb_atoms, 7)-min_energy);
end
%--------------------------------------------------------------------------


function [high_csm_gb_atoms_indx, num_gb_atoms] = find_number_of_gb_atoms(Coords, dump_xyz_lims, CSMcrt)
%Upper Zone
UpperCoords=Coords(Coords(:,4)>=0,1:7);
GB_Coords_Upper=UpperCoords(UpperCoords(:,6)>=CSMcrt,3:5);
Upper_Mean=mean(GB_Coords_Upper(:,2));
%Lower Zone
LowerCoords=Coords(Coords(:,4)<=0,1:7);
GB_Coords_Lower=LowerCoords(LowerCoords(:,6)>=CSMcrt,3:5);
Lower_Mean=mean(GB_Coords_Lower(:,2));

Upper_Split_Line=min((abs(Upper_Mean)+abs(Lower_Mean))/2,0.6*max(Coords(:,4)));
Lower_Split_Line=-Upper_Split_Line;
%Main Zone
Main_GB_Zone_indx = find((Coords(:,4)>=Lower_Split_Line) & (Coords(:,4)<=Upper_Split_Line));
high_csm_gb_atoms_indx = Main_GB_Zone_indx(Coords(Main_GB_Zone_indx, 6)> CSMcrt);
num_gb_atoms = length(high_csm_gb_atoms_indx);
end
%--------------------------------------------------------------------------


function [sorted_vvs_based_on_insertion_prb, insertion_prb] = find_insertion_prb(a, b, CSMcrt, lat_param, lat_type)

if strcmpi(lat_type, 'fcc')
    b_vec_mag = lat_param / sqrt(2);
elseif strcmpi(lat_type, 'bcc')
    b_vec_mag = lat_param * sqrt(3)/2;
end

[Coords_p] = create_box_images(a, b);
% ---------------------------------------------------
Coords = Coords_p(:,3:5);
% ---------------------------------------------------
indx =[];
while length(indx) < 8
    indx_csm = find((Coords_p(:,6) > CSMcrt) & (Coords_p(:,4) < 50) & (Coords_p(:,4) > -50));
    indx = indx_csm(indx_csm <= size(a,1));
    CSMcrt = 0.1 * CSMcrt;
end

if ~isempty(indx)
    success_flag = true;
    % % % %
    r_cut = 3 * b_vec_mag;
    neybs_inds = zeros(size(Coords,1),1);
    cnt3 = 1;
    for cnt2 = 1:length(indx)
        [index,~] = near_neighbor_wo_sort(Coords, indx(cnt2), r_cut);
        if ~isempty(index)
            neybs_inds(cnt3 : cnt3 + length(index)-1) = index;
            cnt3 = cnt3 + length(index);
        end
    end
    neybs_inds(cnt3:end)=[];
    neybs_inds = unique(neybs_inds);
    
    % ----------------------
    unit_w_pad = Coords(neybs_inds, 1:3);
    u_c = indx;
    % ---------------------Triangulation and Voronoi Generation ---------------
    dt = delaunayTriangulation(unit_w_pad);
%     [V, ~] = voronoiDiagram(dt);
%     dt_inner = delaunayTriangulation(Coords(u_c,:));
    % ----
    [CC, r_to_atoms] = circumcenter(dt);
%     cc_inner = tsearchn(Coords(u_c,:), dt_inner.ConnectivityList, CC);
%     indx_cc = find(~isnan(cc_inner));
%     indx_cc = inpolygon();
    [basis_2d_g, base_poly] = basis_from_lammps(b, [0 0 0]);
    [indx_cc1, ~]= inpolygon(CC(:,1),CC(:,3), base_poly(:,1),base_poly(:,3));
    indx_cc2 = CC(:, 2) >= min(Coords(u_c, 2)) & CC(:, 2) <= max(Coords(u_c, 2));
    indx_cc = find(indx_cc1 & indx_cc2);
    
    vv_points = CC(indx_cc, :);
    
    % %%  In-house method for finding the distance between VVs and atoms.
    %     Alternative to Matlab's circumference function.
    %     ====================================================================
    %     tn_inner = tsearchn(Coords(u_c,:), dt_inner.ConnectivityList, V);
    %     indx_v = find(~isnan(tn_inner));
    %     back_map = indx_v;
    %     vv_points = V(indx_v, :);
    %     [voro_to_atom_indx] = voro_to_atom_map(R);
    %     Points=dt.Points;
    %     tmp1 = cellfun(@(x)(x(1))',voro_to_atom_indx,'UniformOutput', false);
    %     close_to_vv_atom_indx = cell2mat(tmp1);
    %     del_vec = V - Points(close_to_vv_atom_indx,:);
    %     dists = sqrt(del_vec(:,1).^2 + del_vec(:,2).^2 + del_vec(:,3).^2);
    % %%  ====================================================================
%     B_vec = lat_param / sqrt(2);
    B_vec = b_vec_mag;
    r = r_to_atoms - B_vec / 2;
    r_of_interest = r(indx_cc);
    [rr, ii] = sort(r_of_interest,'descend');
    insertion_prb = rr ./ sum(rr);
    sorted_vvs_based_on_insertion_prb = vv_points(ii, :);
else
    error('No GBs were found!!!');
end
end
%--------------------------------------------------------------------------


function [Coords_p] = create_box_images(Coords, xyz_lims)
if size(xyz_lims,2)>2
    tilt = xyz_lims(2,3);
    if (tilt<0 & xyz_lims(1,1)<0 )
        tilt_mult = -1;
    else
        tilt_mult = 1;
    end
else
    tilt_mult = 1;
    tilt = 0;
end
% box_size(1:3) = [xyz_lims(1,2)-xyz_lims(1,1)+ tilt_mult * tilt,...
box_size(1:3) = [xyz_lims(1,2)-xyz_lims(1,1)- tilt_mult * tilt,...
    xyz_lims(2,2)-xyz_lims(2,1),...
    xyz_lims(3,2)-xyz_lims(3,1)];


shift_coefs_x =  [0, 1, 1, 0, -1 , -1, -1,  0,  1, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1, 0, 1, 2, 2];
shift_coefs_zt = [0, 0, 1, 1,  1 ,  0, -1, -1, -1, 0, 1, 2, 2, 2, 2, 2, 1, 0,-1,-2,-2,-2,-2,-2,-1];
% shift_coefs_x = [0, 1, 1, 0, -1 , -1, -1, 0, 1];
% shift_coefs_zt = [0, 0, 1, 1, 1 , 0, -1, -1, -1];
Coords_p = repmat(Coords, length(shift_coefs_x), 1);

s1 = size(Coords,1);
s2 = size(Coords,2);
%  creating periodic images around the simulated cell
for i=1:length(shift_coefs_x)
    shift_x = shift_coefs_x(i) * box_size(1) + shift_coefs_zt(i) * (tilt);
    shift_z = shift_coefs_zt(i) * box_size(3);
    if shift_x ~= 0
        Coords_p((i-1)*s1+1:i*s1,3) = Coords(:,3) + shift_x;
    end
    
    if shift_z ~= 0
        Coords_p((i-1)*s1+1:i*s1,5) = Coords(:,5) + shift_z;
    end
end
max_partc_ident =max(Coords(:,1));
Coords_p(s1+1:end, 1) = max_partc_ident+1 : (length(shift_coefs_x)-1) * s1 + max_partc_ident;
end
%--------------------------------------------------------------------------


function [v_to_a_indx] = voro_to_atom_map(R)
all_voro_indx = [R{:}];

%  Building the addresses
addresses  = R;
for i=1:length(addresses)
    addresses{i} = i * ones(size(addresses{i}));
end
addresses_vec = [addresses{:}];

v_to_a_indx = cell(length(all_voro_indx), 1);
for i=1:length(all_voro_indx)
    voro_i_in_R_indx = all_voro_indx == i;
    v_to_a_indx{i} = addresses_vec(voro_i_in_R_indx);
end
v_to_a_indx(find(cell_size(v_to_a_indx,2)==0))=[];
end
% -----------------------------------------------------------------------


function [index,radii] = near_neighbor_wo_sort(input,i,r_cut)
c_unit = input(i,:);
input_centered = input - repmat(c_unit,size(input,1),1);

radii = sqrt(input_centered(:,1).^2 + input_centered(:,2).^2 + input_centered(:,3).^2);

index = find((radii < r_cut) & radii ~= 0);
radii = radii(index);
end
% -----------------------------------------------------------------------


function [index,radii] = near_neighbor_to_coord(input, center_coords, r_cut)
input_centered = input - repmat(center_coords, size(input,1),1);

radii = sqrt(input_centered(:,1).^2 + input_centered(:,2).^2 + input_centered(:,3).^2);
[r_sort,ind_sort] = sort(radii);
index = ind_sort(r_sort<r_cut);
radii = radii(index(2:end));
end
% -----------------------------------------------------------------------


function plot_Bernal_spheres(sph_cent, r_i)
sphere_color = [30 90 200]/255;
[X, Y, Z] = sphere(30);
plot3(sph_cent(1),sph_cent(2),sph_cent(3),'*c','MarkerSize', 20);
h = surf((r_i*X+sph_cent(1)), (r_i*Y+sph_cent(2)), (r_i*Z+sph_cent(3)),...
    'LineStyle', '-', 'FaceColor', sphere_color, 'EdgeColor',...
    'none', 'EdgeAlpha', 0.5, 'LineWidth', 2.0,'FaceAlpha',0.5);
h.SpecularStrength = 0.4;
end
% -----------------------------------------------------------------------


function lammps_data_file_write(tpts_box, box_sizes, out_path)
% ------ removing the atoms just on the edge
tol = 1e-3;
inds_to_chop = (tpts_box(:,4)>=box_sizes(2,2)-tol | tpts_box(:,4)<=box_sizes(2,1)+tol);
tpts_box(inds_to_chop, :) = [];
% ------
fiw = fopen(out_path,'w');
cnt = 1;
line{cnt} = ['#LAMMPS data file generated by Arash Banadaki\n'];cnt = cnt + 1;
line{cnt} = [num2str(size(tpts_box,1)),' atoms\n'];cnt = cnt + 1;
line{cnt} = ['2   atom types\n'];cnt = cnt + 1;
line{cnt} = ['\n'];cnt = cnt + 1;
if size(box_sizes, 2) == 2
    line{cnt} = [num2str(box_sizes(1,1),'%14.10e'),' ',num2str(box_sizes(1,2),'%14.10e'),' xlo xhi\n'];cnt = cnt + 1;
    line{cnt} = [num2str(box_sizes(2,1),'%14.10e'),' ',num2str(box_sizes(2,2),'%14.10e'),' ylo yhi\n'];cnt = cnt + 1;
    line{cnt} = [num2str(box_sizes(3,1),'%14.10e'),' ',num2str(box_sizes(3,2),'%14.10e'),' zlo zhi\n'];cnt = cnt + 1;
elseif size(box_sizes, 2) == 3
    if box_sizes(1,1)<0
        if box_sizes(1,1) == box_sizes(2,3)
            line{cnt} = [num2str(0),' ',num2str(box_sizes(1,2),'%14.10e'),' xlo xhi\n'];cnt = cnt + 1;
        else
            line{cnt} = [num2str(box_sizes(1,1),'%14.10e'),' ',num2str(box_sizes(1,2)-abs(box_sizes(2,3)),'%14.10e') ,' xlo xhi\n'];cnt = cnt + 1;
        end
    else
        line{cnt} = [num2str(box_sizes(1,1),'%14.10e'),' ',num2str(box_sizes(1,2)-box_sizes(2,3),'%14.10e') ,' xlo xhi\n'];cnt = cnt + 1;
        %         line{cnt} = [num2str(box_sizes(2,3)),' ',num2str(box_sizes(1,2)) ,' xlo xhi\n'];cnt = cnt + 1;
    end
    line{cnt} = [num2str(box_sizes(2,1),'%14.10e'),' ',num2str(box_sizes(2,2),'%14.10e'),' ylo yhi\n'];cnt = cnt + 1;
    line{cnt} = [num2str(box_sizes(3,1),'%14.10e'),' ',num2str(box_sizes(3,2),'%14.10e'),' zlo zhi\n'];cnt = cnt + 1;
    line{cnt} = [num2str(box_sizes(1,3),'%14.10e'),' ',num2str(box_sizes(2,3),'%14.10e'),' ',num2str(box_sizes(3,3),'%14.10e'),' xy xz yz\n'];cnt = cnt + 1;
else
    error('box sizes can not be identified!');
end
%     line{cnt} = ['\n'];cnt = cnt + 1;
%     line{cnt} = ['Masses\n'];cnt = cnt + 1;
%     line{cnt} = ['\n'];cnt = cnt + 1;
%     line{cnt} = ['1	26.98\n'];cnt = cnt + 1;
line{cnt} = ['\n'];cnt = cnt + 1;
line{cnt} = ['Atoms\n'];cnt = cnt + 1;
line{cnt} = ['\n'];
for i = 1:cnt
    fprintf(fiw, line{i});
end

mat_print = tpts_box;
fprintf(fiw, '%g\t%g\t%g\t%g\t%g\n', reshape(mat_print', 5, size(tpts_box, 1)));
fclose(fiw);
end
%--------------------------------------------------------------------------


function [] = add_single_crystal_block(elem, lat_param, csm_crit, gb_main_folder, l_p_go_attr, cnt1)

if strcmpi(elem,'al')
    min_energy = -3.36000010077;
    lat_param = 4.05;
    lat_type = 'fcc';
    l_p_po = lat_param * [[0 0.5 0.5]', [0.5 0 0.5]', [0.5 0.5 0]'];
elseif strcmpi(elem,'cu')
    min_energy = -3.53999996838;
    lat_param = 3.615;
    lat_type = 'fcc';
    l_p_po = lat_param * [[0 0.5 0.5]', [0.5 0 0.5]', [0.5 0.5 0]'];
elseif strcmpi(elem,'ni')
    min_energy = -4.45000000527;
    lat_param = 3.52;
    lat_type = 'fcc';
    l_p_po = lat_param * [[0 0.5 0.5]', [0.5 0 0.5]', [0.5 0.5 0]'];
elseif strcmpi(elem,'fe')
    min_energy = -4.1224310073005;
    lat_param = 2.855324;    
    lat_type = 'bcc';
    l_p_po = lat_param * [[-0.5 0.5 0.5]', [0.5 -0.5 0.5]', [0.5 0.5 -0.5]'];
elseif strcmpi(elem,'mo')
    min_energy = -6.8200023329166;
    lat_param = 3.1472;    
    lat_type = 'bcc';
    l_p_po = lat_param * [[-0.5 0.5 0.5]', [0.5 -0.5 0.5]', [0.5 0.5 -0.5]'];
end


% gb_dump_file = fullfile(gb_main_folder, ['dump_' num2str(cnt1)]);%%%% origianl, must come back and check the next line
gb_dump_file = fullfile(gb_main_folder, ['dump_' num2str(cnt1+1)]);
bpn_po = bpn_extractor(gb_dump_file);
bpn_p = (l_p_po/lat_param)^-1 * bpn_po;
int_pln_spc = int_planar_spc(bpn_p, l_p_po);
int_spc = floor(lat_param / int_pln_spc) * int_pln_spc;
% int_spc = int_pln_spc;
% ------
[Coords, dump_xyz_lims] = f_read(gb_dump_file);
[gb_attr] = find_gb_attributes(Coords, csm_crit);

if gb_attr.mean_gb_y >= 0
    y_axis_extent_go = [0 dump_xyz_lims(2, 2) 0]';
    cent_single_cryst_region = mean(Coords(gb_attr.upper_single_cryst_index, 3:5));
    l_p_go = l_p_go_attr.upper;
else
    y_axis_extent_go = [0 dump_xyz_lims(2, 1) 0]';
    cent_single_cryst_region = mean(Coords(gb_attr.lower_single_cryst_index, 3:5));
    l_p_go = l_p_go_attr.lower;
end

% figure; plot_gb(Coords, gb_attr, int_spc);view([1 0 0]);
% --------------------------
% g_basis = [[1 0 0]' [0 1 13]' [0 -13 1]'];
% g_basis_n = [g_basis(:,1)/norm(g_basis(:,1)),...
%              g_basis(:,2)/norm(g_basis(:,2)),...
%              g_basis(:,3)/norm(g_basis(:,3))]';

% % % l_po_go  = (l_p_po/lat_param) * l_p_go^-1;
l_po_go  = (l_p_po * l_p_go^-1)^-1;

% 
% l_po_go = [l_po_go(:,1)/norm(l_po_go(:,1)),...
%     l_po_go(:,2)/norm(l_po_go(:,2)),...
%     l_po_go(:,3)/norm(l_po_go(:,3))];
% g_basis_n = l_po_go';
g_basis_n = l_po_go^-1;
% --------------------------
% y_axis_extent_po = g_basis_n' * y_axis_extent_go;
y_axis_extent_po = g_basis_n * y_axis_extent_go;
num_supercells_in_po = ceil(max(abs(y_axis_extent_po))/lat_param) * 2;
x_num_cells = num_supercells_in_po;   y_num_cells = num_supercells_in_po;  z_num_cells = num_supercells_in_po; % number of cells to be generated initially
c1 = gen_atom(lat_param, lat_type, x_num_cells, y_num_cells, z_num_cells); % Generatin the FCC structure
% C1 = (g_basis_n * c1')';
C1 = (l_po_go * c1')';
% ------
[index, ~] = near_neighbor_to_coord(c1, mean(c1), 10);
cent_indx = index(1);
cent_coord = C1(cent_indx, :);

[index1, ~] = near_neighbor_to_coord(Coords(:,3:5), cent_single_cryst_region, 10);
cent_sing_cryst_atom = Coords(index1(1), 3:5);
C1 = C1 - repmat(cent_coord-cent_sing_cryst_atom, size(C1,1), 1);
% ------
[basis_2d_g, base_poly] = basis_from_lammps(dump_xyz_lims, [0 0 0]);
[in_p on_p]= inpolygon(C1(:,1), C1(:,3), base_poly(:,1), base_poly(:,3));
% ------
tmp = C1(in_p, :);
if gb_attr.mean_gb_y >= 0
    C2 = tmp(tmp(:,2) >= cent_sing_cryst_atom(2), :);
    C3 = Coords(Coords(:,4) <= cent_sing_cryst_atom(2), 3:5);
else
    C2 = tmp(tmp(:,2) < cent_sing_cryst_atom(2), :);
    C3 = Coords(Coords(:,4) >= cent_sing_cryst_atom(2), 3:5);
end
C = [C2;C3];
C(:,2) = C(:,2) - gb_attr.mean_gb_y;
to_chop_in_y_direction = C(:, 2) >= dump_xyz_lims(2,2) | C(:, 2) <= dump_xyz_lims(2,1);
C(to_chop_in_y_direction, :) = [];

[~, unq_rows] = UniqueRows_Tol(C(:, 1:3), 1e-1 * lat_param);
C = C(unq_rows, :);

% plot3(C2(:, 1), C2(:, 2), C2(:, 3),'^k');axis equal; axis off;
% plot3(C3(:, 1), C3(:, 2), C3(:, 3),'^k');axis equal; axis off;
% figure; plot3(C(:, 1), C(:, 2), C(:, 3),'ok');axis equal; axis off;view([1 0 0]);
% ------
lammps_data(:, 1) = 1:size(C, 1);
upper_indx = (C(:, 2) >= 0);
lower_indx = (C(:, 2) < 0);
lammps_data(upper_indx, 2) = 1;
lammps_data(lower_indx, 2) = 2;
lammps_data(:, 3:5) = C;
% % ------
lammps_data_file_name  = fullfile(gb_main_folder, ['dump_' num2str(cnt1+1) '.t']);
lammps_data_file_write(lammps_data, dump_xyz_lims, lammps_data_file_name);
gb_eng_path = fullfile(gb_main_folder,'GB_Eng.dat');
% ----------------------------------
next_dump = fullfile(gb_main_folder, ['dump_' num2str(cnt1+1)]);
pot_file = ['pot_', elem];
%cmnd1= 'mpirun -np 10 lammps-daily -log none -var';
global lammps_exe_path;
cmnd1= ['mpirun ', lammps_exe_path, ' -log none -var'];
cmnd2= [' potname1 ', pot_file,' -var potname2 ', pot_file,' -var '];
%         cmnd3 = ['fname1 ',lammps_data_file_name,' -var fname2 ',copper_dump, ' < in.GB_convert_AtoB.txt'];
cmnd3 = ['fname1 ',lammps_data_file_name,' -var fname2 ',next_dump,...
    ' -var fname3 ',gb_eng_path, ' -var min_eng ', num2str(min_energy, 15),...
    ' -var lat_type ', lat_type, ' < in.minimize'];
cmnd = [cmnd1 cmnd2 cmnd3];
system(cmnd);
end
%--------------------------------------------------------------------------


function [basis_2d_g, base_poly] = basis_from_lammps(b, origin)
if size(b,2)==3
    t1 = [abs(b(1,2)-b(1,1))-abs(b(2,3)), 0, 0];
    t2 = [b(2,3), 0, b(3,2)-b(3,1)];
    basis_2d_g = [t1', t2'];
    base_poly(1, :) = origin;
    base_poly(2, :) = base_poly(1, :) + t1;
    base_poly(3, :) = base_poly(2, :) + t2;
    base_poly(4, :) = base_poly(1, :) + t2;
else
    t1 = [abs(b(1,2)-b(1,1)), 0, 0];
    t2 = [0, 0, b(3,2)-b(3,1)];
    basis_2d_g = [t1', t2'];
    base_poly(1, :) = origin;
    base_poly(2, :) = base_poly(1, :) + t1;
    base_poly(3, :) = base_poly(2, :) + t2;
    base_poly(4, :) = base_poly(1, :) + t2;
end
base_poly = [base_poly;base_poly(1,:)];
end
%--------------------------------------------------------------------------


function S = gen_atom(lat_param, lat_type, x_size,y_size,z_size)
if strcmpi(lat_type, 'fcc')
    b1 = [0,0,0]; b2 = [0,1,1]/2; b3 = [1,0,1]/2; b4 = [1,1,0]/2;
    unit_cell = [b1;b2;b3;b4];
    % ---Repeating in global X direction
    S = repmat(unit_cell,x_size+1,1);
    delta_x = reshape(repmat(0:x_size,4,1),1,4*(x_size+1));
    S(:,1) = S(:,1)+ delta_x';
    % ---Repeating in global Y direction
    delta_y = reshape(repmat(0:y_size,size(S,1),1),1,size(S,1)*(y_size+1));
    S = repmat(S,y_size+1,1);
    S(:,2) = S(:,2)+ delta_y';
    % ---Repeating in global Z direction
    delta_z = reshape(repmat(0:z_size,size(S,1),1),1,size(S,1)*(z_size+1));
    S = repmat(S,z_size+1,1);
    S(:,3) = S(:,3)+ delta_z';
    % -----------------
    S(S(:,1)>x_size,:)=[];
    S(S(:,2)>y_size,:)=[];
    S(S(:,3)>z_size,:)=[];
    S = S * lat_param;
    S(:,1)=S(:,1)-x_size/2; S(:,2)=S(:,2)-y_size/2; S(:,3)=S(:,3)-z_size/2; % Moving the center to the origin
elseif strcmpi(lat_type, 'bcc')
    b1 = [0,0,0]; b2 = [1,1,1]/2;
    unit_cell = [b1;b2];
    % ---Repeating in global X direction
    S = repmat(unit_cell,x_size+1,1);
    delta_x = reshape(repmat(0:x_size,2,1),1,2*(x_size+1));
    S(:,1) = S(:,1)+ delta_x';
    % ---Repeating in global Y direction
    delta_y = reshape(repmat(0:y_size,size(S,1),1),1,size(S,1)*(y_size+1));
    S = repmat(S,y_size+1,1);
    S(:,2) = S(:,2)+ delta_y';
    % ---Repeating in global Z direction
    delta_z = reshape(repmat(0:z_size,size(S,1),1),1,size(S,1)*(z_size+1));
    S = repmat(S,z_size+1,1);
    S(:,3) = S(:,3)+ delta_z';
    % -----------------
    S(S(:,1)>x_size,:)=[];
    S(S(:,2)>y_size,:)=[];
    S(S(:,3)>z_size,:)=[];
    S = S * lat_param;
    S(:,1)=S(:,1)-x_size/2; S(:,2)=S(:,2)-y_size/2; S(:,3)=S(:,3)-z_size/2; % Moving the center to the origin
end
end
%--------------------------------------------------------------------------


function plot_gb(Coords_i, gb_attr_i, int_pln_spc)
plot3(Coords_i(:,3),Coords_i(:,4),Coords_i(:,5),'.b');axis equal;axis off;hold on;
plot3(Coords_i(gb_attr_i.indx,3),Coords_i(gb_attr_i.indx,4),Coords_i(gb_attr_i.indx,5),'or');
% l_i = find(Coords_i(:,4) > gb_attr_i.bounds(2)+1*int_pln_spc-1e-3);
l_i = find(Coords_i(:,4) > gb_attr_i.bounds(2));
% plot3(Coords_i(l_i,3),Coords_i(l_i,4),Coords_i(l_i,5),'og');
end
%--------------------------------------------------------------------------


function [pln_spc] = int_planar_spc(bpn_p, l_p_po)
l_rp_po = reciprocal_mat(l_p_po);
l_po_rp = l_rp_po ^ -1;
l_p_rp = l_po_rp * l_p_po;
bpn_rp = l_p_rp * bpn_p;

bpn_rp = RationalFinder(bpn_rp);
pln_spc = 1/norm(l_rp_po * bpn_rp);
end
%--------------------------------------------------------------------------


function [rl_g_go] = reciprocal_mat(l_g_go)
L3 = cross(l_g_go(:, 1), l_g_go(:, 2)') / det(l_g_go);
L1 = cross(l_g_go(:, 2), l_g_go(:, 3)') / det(l_g_go);
L2 = cross(l_g_go(:, 3), l_g_go(:, 1)') / det(l_g_go);
rl_g_go = [L1; L2; L3]';
end
%--------------------------------------------------------------------------


function [bpn] = bpn_extractor(folder_name)
tmp = strsplit(folder_name,'/');
gb_name = tmp{end-1};
n1 = strfind(gb_name,'N1');
n2 = strfind(gb_name,'N2');
tmp = strsplit(gb_name(n1:n2) ,'_');
bpn(1,1) = str2num(tmp{2});
bpn(2,1) = str2num(tmp{3});
bpn(3,1) = str2num(tmp{4});
end
%--------------------------------------------------------------------------


function [gb_attributes] = find_gb_attributes(Coords, CSMcrt)
%Upper Zone
UpperCoords=Coords(Coords(:,4)>=0,1:7);
GB_Coords_Upper=UpperCoords(UpperCoords(:,6)>=CSMcrt,3:5);
Upper_Mean=mean(GB_Coords_Upper(:,2));
%Lower Zone
LowerCoords=Coords(Coords(:,4)<=0,1:7);
GB_Coords_Lower=LowerCoords(LowerCoords(:,6)>=CSMcrt,3:5);
Lower_Mean=mean(GB_Coords_Lower(:,2));

Upper_Split_Line=min((abs(Upper_Mean)+abs(Lower_Mean))/2,0.6*max(Coords(:,4)));
Lower_Split_Line=-Upper_Split_Line;
%Main Zone
Main_GB_Zone_indx = find((Coords(:,4)>=Lower_Split_Line) & (Coords(:,4)<=Upper_Split_Line));
high_csm_gb_atoms_indx = Main_GB_Zone_indx(Coords(Main_GB_Zone_indx, 6)> CSMcrt);
num_gb_atoms = length(high_csm_gb_atoms_indx);
% ----
mean_gb_depth = sum(Coords(high_csm_gb_atoms_indx, 4) .* ...
    Coords(high_csm_gb_atoms_indx, 6)) / sum(Coords(high_csm_gb_atoms_indx, 6));
[gb_ymax, gb_ymax_indx]= max(Coords(high_csm_gb_atoms_indx, 4));
[gb_ymin, gb_ymin_indx]= min(Coords(high_csm_gb_atoms_indx, 4));
gb_bounds = [gb_ymin gb_ymax];
gb_bounds_indx = [high_csm_gb_atoms_indx(gb_ymin_indx) high_csm_gb_atoms_indx(gb_ymax_indx)];
% ----
upper_single_cryst_region_index = find(Coords(:,4)> mean_gb_depth & Coords(:, 6)< CSMcrt);
lower_single_cryst_region_index = find(Coords(:,4)< mean_gb_depth & Coords(:, 6)< CSMcrt);
% ----
gb_attributes.bounds = gb_bounds;
gb_attributes.bounds_indx = gb_bounds_indx;
gb_attributes.num_atoms = num_gb_atoms;
gb_attributes.indx = high_csm_gb_atoms_indx;
gb_attributes.mean_gb_y = mean_gb_depth;
gb_attributes.lower_single_cryst_index  = lower_single_cryst_region_index;
gb_attributes.upper_single_cryst_index = upper_single_cryst_region_index;
end
%--------------------------------------------------------------------------

function [l_p_go] = l_p_go_finder(gb_data, gb_attr, csm_crit, lat_param, lat_type)
% Coords_in = gb_data.Coords;
% dump_xyz_lims = gb_data.dump_xyz_lims;
% 
% if gb_attr.mean_gb_y >= 0
%     upper_lower = 1;
% else
%     upper_lower = -1;
% end
% 
% [Coords] = create_box_images(Coords_in, dump_xyz_lims);
% [single_crystal_lower_y, single_crystal_upper_y] = single_crystal_zone(Coords, csm_crit);
% if upper_lower == 1
%     single_crystal_mean_y = single_crystal_upper_y;
% elseif upper_lower == -1
%     single_crystal_mean_y = single_crystal_lower_y;
% end
% plt_flag = false;
% l_p_go = find_first_atom(Coords, single_crystal_mean_y, lat_param, plt_flag);
Coords_in = gb_data.Coords;
dump_xyz_lims = gb_data.dump_xyz_lims;

[Coords] = create_box_images(Coords_in, dump_xyz_lims);
[single_crystal_lower_y, single_crystal_upper_y] = single_crystal_zone(Coords, csm_crit);

plt_flag = false;
l_p_go_upper = find_first_atom(Coords, single_crystal_upper_y, lat_param, plt_flag, lat_type);
l_p_go_lower = find_first_atom(Coords, single_crystal_lower_y, lat_param, plt_flag, lat_type);
l_p_go.lower = l_p_go_lower;
l_p_go.upper = l_p_go_upper;
end
% -------------------------------------------------------------------------


% function [Lower_Mean, Upper_Mean] = single_crystal_zone (Coords, CSMcrt)
function [Lower_Split_Line, Upper_Split_Line] = single_crystal_zone(Coords, CSMcrt)
%Upper Zone
CSMcrt = 0.01;
UpperCoords = Coords(Coords(:,4)>=0,1:7);
GB_Coords_Upper = UpperCoords(UpperCoords(:,6)>=CSMcrt,3:5);
Upper_Mean = mean(GB_Coords_Upper(:,2));
%Lower Zone
LowerCoords=Coords(Coords(:,4)<=0,1:7);
GB_Coords_Lower=LowerCoords(LowerCoords(:,6)>=CSMcrt,3:5);
Lower_Mean = mean(GB_Coords_Lower(:,2));
% ---------
Upper_Split_Line=min((abs(Upper_Mean)+abs(Lower_Mean))/2,0.6*max(Coords(:,4)));
Lower_Split_Line = -Upper_Split_Line;
end
% -------------------------------------------------------------------------


function [l_p_go] = find_first_atom(Coords, single_crystal_mean_y, lat_param, plt_flag, lat_type)

estimated_first_atom_coords = [mean(Coords(:,3)), single_crystal_mean_y,...
    mean(Coords(:,5))];
vec = Coords(:,3:5) - repmat(estimated_first_atom_coords, size(Coords,1), 1);
[~, first_atom_indx] = min(sqrt(vec(:,1).^2 + vec(:,2).^2 + vec(:,3).^2));

if strcmpi(lat_type, 'fcc')
    burgers = lat_param / sqrt(2);
    r_cut = 1.005 * burgers;
    combos = nchoosek(1:12,3);
elseif strcmpi(lat_type, 'bcc')
    burgers = lat_param * sqrt(3)/ 2;
    r_cut = 1.005 * burgers;
    combos = nchoosek(1:8,3);
end

[indx_2, radii_2] = near_neighbor(Coords(:,3:5), first_atom_indx, r_cut);
all_neighs_indx = indx_2(2:end);
possible_triangles_indx = all_neighs_indx(combos);

first_vec = Coords(possible_triangles_indx(:,1),3:5)-Coords(possible_triangles_indx(:,2),3:5);
first_vec_mag = sqrt(first_vec(:,1).^2 + first_vec(:,2).^2 + first_vec(:,3).^2);

second_vec = Coords(possible_triangles_indx(:,1),3:5)-Coords(possible_triangles_indx(:,3),3:5);
second_vec_mag = sqrt(second_vec(:,1).^2 + second_vec(:,2).^2 + second_vec(:,3).^2);

third_vec = Coords(possible_triangles_indx(:,2),3:5)-Coords(possible_triangles_indx(:,3),3:5);
third_vec_mag = sqrt(third_vec(:,1).^2 + third_vec(:,2).^2 + third_vec(:,3).^2);

if strcmpi(lat_type, 'fcc')
    dist_diff = abs([first_vec_mag second_vec_mag third_vec_mag]) - radii_2(1);
    mag_dist_diff = sqrt(dist_diff(:,1).^2 + dist_diff(:,2).^2 +dist_diff(:,3).^2);
elseif strcmpi(lat_type, 'bcc')
	dist_diff = [first_vec_mag second_vec_mag third_vec_mag];
    mag_dist_diff = abs(dist_diff(:,1) + dist_diff(:,2) + dist_diff(:,3)-(2 * lat_param + sqrt(2)*lat_param));
end
% equilat_tri_indx = find(mag_dist_diff < 1e-3);
[min_dist ia] = sort(mag_dist_diff);
if min_dist(1) < 1e-2 % This precision is too low. However the current dump files have
    equilat_tri_indx = ia(1);
else
    error('No basis has been found. The found atoms don''t belong to the single crystal!');
end
vertices  = possible_triangles_indx(equilat_tri_indx(1),:);

% % % %
if length(vertices)==3
    basis_vec_1 = (Coords(vertices(1),3:5)-Coords(first_atom_indx,3:5))';
    %     basis_vec_1 = basis_vec_1 / norm(basis_vec_1);
    basis_vec_2 = (Coords(vertices(2),3:5)-Coords(first_atom_indx,3:5))';
    %     basis_vec_2 = basis_vec_2 / norm(basis_vec_2);
    basis_vec_3 = (Coords(vertices(3),3:5)-Coords(first_atom_indx,3:5))';
    %     basis_vec_3 = basis_vec_3 / norm(basis_vec_3);
    if det([basis_vec_1 basis_vec_2 basis_vec_3]) > 0
        l_p_go = [basis_vec_1 basis_vec_2 basis_vec_3];
    else
        l_p_go = [-basis_vec_1 -basis_vec_2 -basis_vec_3];
    end
else
    error('The closest nn set is not consistent.');
end
%
if plt_flag
    plot3(Coords(all_neighs_indx,3),Coords(all_neighs_indx,4),Coords(all_neighs_indx,5),'ob');axis equal;hold on;
    plot3(Coords(first_atom_indx,3),Coords(first_atom_indx,4),Coords(first_atom_indx,5),'*r');hold on ;axis equal;
    plot3(Coords(vertices,3),Coords(vertices,4),Coords(vertices,5),'*r');axis equal;
    quiver3(repmat(Coords(first_atom_indx,3),1,3),repmat(Coords(first_atom_indx,4),1,3),repmat(Coords(first_atom_indx,5),1,3), l_p_go(1,:),l_p_go(2,:),l_p_go(3,:),0)
    % hold off;
end
end
% -----------------------------------------------------------------------


function [index,radii] = near_neighbor(input,i,r_cut)
c_unit = input(i,:);
input_centered = input - repmat(c_unit,size(input,1),1);
% radii = cell2mat(arrayfun(@(x) norm(input_centered(x,:)),1:size(input_centered,1),'uni',false));
radii = sqrt(input_centered(:,1).^2 + input_centered(:,2).^2 + input_centered(:,3).^2);
[r_sort,ind_sort] = sort(radii);
index = ind_sort(r_sort<r_cut);
radii = radii(index(2:end));
end
% -----------------------------------------------------------------------


function [] = in_plane_shift(MC_data, gb_main_folder, csm_crit, cnt1)

elem = MC_data.elem;
if strcmpi(elem,'al')
    min_energy = -3.36000010077;
    lat_param = 4.05;
    lat_type = 'fcc';
elseif strcmpi(elem,'cu')
    min_energy = -3.53999996838;
    lat_param = 3.615;
    lat_type = 'fcc';
elseif strcmpi(elem,'ni')
    min_energy = -4.45000000527;
    lat_param = 3.52;
    lat_type = 'fcc';
elseif strcmpi(elem,'fe')
    min_energy = -4.1224310073005;
    lat_param = 2.855324;    
    lat_type = 'bcc';
elseif strcmpi(elem,'mo')
    min_energy = -6.8200023329166;
    lat_param = 3.1472;    
    lat_type = 'bcc';
end

Coords = MC_data.Coords;
dump_xyz_lims = MC_data.dump_xyz_lims;
[gb_attr] = find_gb_attributes(Coords, csm_crit);

% % ------
up_cryst_ind = find(Coords(:,4) > gb_attr.mean_gb_y);
[basis_2d_g, base_poly] = basis_from_lammps(dump_xyz_lims, [0 0 0]);
x_transl = rand() * basis_2d_g(:, 1)';
z_transl = rand() * basis_2d_g(:, 2)';
if strcmpi(lat_type, 'fcc')
    y_transl = lat_param/sqrt(2) * [0 1 0];
elseif strcmpi(lat_type, 'bcc')
    y_transl = lat_param *sqrt(3)/2 * [0 1 0];
end
transl_vec = x_transl + y_transl + z_transl;
Coords(up_cryst_ind, 3:5) = Coords(up_cryst_ind, 3:5) + repmat(transl_vec, length(up_cryst_ind), 1);
dump_xyz_lims(2,2) = dump_xyz_lims(2,2) + y_transl(2);
% % ------
% % wrapping the periodic BCs
[Coords_im] = create_box_images(Coords, dump_xyz_lims);
[in_p on_p]= inpolygon(Coords_im(:,3), Coords_im(:,5), base_poly(:,1), base_poly(:,3));
% Coords_im =Coords;
% in_p = ones(size(Coords,1),1);
% % ------
lammps_data(:, 1) = 1:length(find(in_p));
upper_indx = (Coords_im(in_p, 4) >= 0);
lower_indx = (Coords_im(in_p, 4) < 0);
lammps_data(upper_indx, 2) = 1;
lammps_data(lower_indx, 2) = 2;
lammps_data(:, 3:5) = Coords_im(in_p, 3:5);
% % ------
dump_file = fullfile(gb_main_folder, ['dump_' num2str(cnt1)]);
lammps_data_file_name = [dump_file, '.t'];
gb_eng_path = fullfile(gb_main_folder,'GB_Eng.dat');
% lammps_data_file_write(Coords_im(in_p, 1:5), dump_xyz_lims, [dump_file '.t']);
lammps_data_file_write(lammps_data, dump_xyz_lims, [dump_file '.t']);
% ----------------------------------
length_cnt1 = length(num2str(cnt1));
next_dump = [dump_file(1:end-length_cnt1) num2str(cnt1+1)];
pot_file = ['pot_', elem];
%cmnd1= 'mpirun -np 10 lammps-daily -log none -var';
global lammps_exe_path;
cmnd1= ['mpirun ', lammps_exe_path, ' -log none -var'];
cmnd2= [' potname1 ', pot_file,' -var potname2 ', pot_file,' -var '];
%         cmnd3 = ['fname1 ',lammps_data_file_name,' -var fname2 ',copper_dump, ' < in.GB_convert_AtoB.txt'];
cmnd3 = ['fname1 ',lammps_data_file_name,' -var fname2 ',next_dump,...
         ' -var fname3 ',gb_eng_path, ' -var min_eng ', num2str(min_energy, 15),...
         ' -var lat_type ', lat_type, ' < in.minimize'];
cmnd = [cmnd1 cmnd2 cmnd3];
system(cmnd);
end

% -----------------------------------------------------------------------


function [] = move_back(lat_param, gb_main_folder, cnt1)
dump_file = fullfile(gb_main_folder, ['dump_' num2str(cnt1+1)]);
[Coords, dump_xyz_lims] = f_read(dump_file);
mid_sim_box = mean([min(Coords(:,4)), max(Coords(:,4))]);
Coords(:,4) = Coords(:,4) - mid_sim_box;
dump_xyz_lims(2,1:2) = [min(Coords(:,4))-1e-3, max(Coords(:,4))+1e-3];
% ------
fiw = fopen(dump_file,'w');
cnt = 1;
line{cnt} = ['ITEM: TIMESTEP\n'];cnt = cnt + 1;
line{cnt} = ['0\n'];cnt = cnt + 1;
line{cnt} = ['ITEM: NUMBER OF ATOMS\n'];cnt = cnt + 1;
line{cnt} = [num2str(size(Coords,1)),'\n'];cnt = cnt + 1;
if size(dump_xyz_lims,2)>2
    line{cnt} = ['ITEM: BOX BOUNDS xy xz yz pp ff pp\n'];cnt = cnt + 1;
    line{cnt} = [num2str(dump_xyz_lims(1,1)),' ',num2str(dump_xyz_lims(1,2)),' ',num2str(dump_xyz_lims(1,3)),'\n'];cnt = cnt + 1;
    line{cnt} = [num2str(dump_xyz_lims(2,1)),' ',num2str(dump_xyz_lims(2,2)),' ',num2str(dump_xyz_lims(2,3)),'\n'];cnt = cnt + 1;
    line{cnt} = [num2str(dump_xyz_lims(3,1)),' ',num2str(dump_xyz_lims(3,2)),' ',num2str(dump_xyz_lims(3,3)),'\n'];cnt = cnt + 1;
else
    line{cnt} = ['ITEM: BOX BOUNDS pp ff pp\n'];cnt = cnt + 1;
    line{cnt} = [num2str(dump_xyz_lims(1,1)),' ',num2str(dump_xyz_lims(1,2)),'\n'];cnt = cnt + 1;
    line{cnt} = [num2str(dump_xyz_lims(2,1)),' ',num2str(dump_xyz_lims(2,2)),'\n'];cnt = cnt + 1;
    line{cnt} = [num2str(dump_xyz_lims(3,1)),' ',num2str(dump_xyz_lims(3,2)),'\n'];cnt = cnt + 1;
end
line{cnt} = ['ITEM: ATOMS id type x y z c_csym c_eng','\n'];
for i = 1:cnt
    fprintf(fiw, line{i});
end
% mat_print = [(1:size(tpts_box,1))', ones(size(tpts_box, 1), 1), tpts_box];
mat_print = Coords;
fprintf(fiw, '%g\t%g\t%g\t%g\t%g\t%g\t%g\n', reshape(mat_print', 7, size(Coords, 1)));
fclose(fiw);
end
% -----------------------------------------------------------------------


function [gb_area_uc] = reag_gb_area(file_name)
fid = fopen(file_name, 'r');
line = fgetl(fid); fclose(fid);
gb_area_uc = str2num(line);
end
