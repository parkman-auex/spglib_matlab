%% note: POSCAR is Direct mode
% input : Stucture information
%                             a_crystal_constance
%                             a1,a2,a3,Rm=[a1;a2;a3]
%                             atom_name;atom_num
%                             sites
%          element_information:
%                             sites 
% Output : POSCAR: 
%                             a_crystal_constance
%                             a1,a2,a3,Rm=[a1;a2;a3]
%                             atom_name;atom_num
%                             coordinates pattern

% note : site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
% Rm
    %a1 a1x a1y a1z
    %a2 a2x a2y a2z
    %a3 a3x a3y a3z
% unit A?e-6
% usage : POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename);






function [Rm,sites] = POSCAR_gen(Rm,sites,Atom_name,Atom_num,filename)
%% nargin
    if nargin < 5
        filename = "POSCAR_gen";
    end
%% init struct
    %site=struct('seq',[],'inseq',[],'rc1',[],'rc2',[],'rc3',[],'name',[],'nameseq',[]);
    title = "POSCAR Gen by MATLAB";
%% write POSCAR
% Initialize variables
    
    %filename = "POSCAR_"+num2str(term2)+".vasp";
    fileID = fopen(filename,'w');
%% Rm
    a_crystal_constance=1;
%% 
    fprintf(fileID,"%s\n",title);
    fprintf(fileID,"%d\n",a_crystal_constance);
    %fprintf(fileID,"  ",Rm(i,j));
    for i=1:3
        for j=1:3
                fprintf(fileID,"  %f",Rm(i,j));
        end
        fprintf(fileID,"\n");
    end
    for i=1:length(Atom_name)
        fprintf(fileID,"%s ",Atom_name(i));
    end
    fprintf(fileID,"\n  ");
    for i=1:length(Atom_num)
        fprintf(fileID,"%d ",Atom_num(i));
    end
    fprintf(fileID,"\n");
    fprintf(fileID,"Direct\n  ");
% sites
    [~,sites_num]=size(sites);
    for i=1:sites_num
        fprintf(fileID,"%f  ",regular_rc(sites(i).rc1));
        fprintf(fileID,"%f  ",regular_rc(sites(i).rc2));
        fprintf(fileID,"%f  ",regular_rc(sites(i).rc3));
%         if ~strcmp(string(sites(i).name),"")
%             %fprintf(fileID,"%s\n  ",sites(i).name);
%             fprintf(fileID,"\n  ");
%         else
%             fprintf(fileID,"\n  ");
%         end
        fprintf(fileID,"\n  ");
    end
end
function rc = regular_rc(Rc)
    rc = mod(Rc,1);
end
