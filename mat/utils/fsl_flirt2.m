function [ out_image, R_mat] = fsl_flirt2(path,ref_image,in_inmage,DOF,res,init,apply)
%fsl fast function
%   filename : string filename WITHOUT nii, f: f value for fast
switch nargin
    case 2
        DOF = 12;
        res = [1,1,2];
        init = false;
        apply = false;
    case 3
        res = [1, 1, 2];
        init = false;
        apply = false;
    case 4
        init = false;
        apply = false;
    case 5
        init = false;
        apply = false;
        
end
[~, ~, ~] = mkdir([path 'forFSL']);
save_nii(make_nii((rot90(in_inmage,-1)),res),[path 'forFSL/temp1.nii'])
save_nii(make_nii((rot90(ref_image,-1)),res),[path 'forFSL/temp2.nii'])
if init
    if apply
        command = ['flirt -in ' path 'forFSL/temp1 -ref ' path 'forFSL/temp2 -dof ' num2str(DOF) ' -applyxfm -out ' path 'forFSL/temp3  -init ' path 'forFSL/R_mat -interp sinc'];
    else
        command = ['flirt -in ' path 'forFSL/temp1 -ref ' path 'forFSL/temp2 -dof ' num2str(DOF) ' -out ' path 'forFSL/temp3  -init ' path 'forFSL/R_mat -interp sinc'];
    end
else
    command = ['flirt -in ' path 'forFSL/temp1 -ref ' path 'forFSL/temp2 -out ' path 'forFSL/temp3 -omat ' path 'forFSL/R_mat -dof ' num2str(DOF)];
end
system(command);
load([path 'forFSL/R_mat'],'-ascii','R_mat');

gunzip([path 'forFSL/*.nii.gz']);
imstr = load_nii([path 'forFSL/temp3.nii']);
out_image = rot90(imstr.img);

% delete([path 'forFSL/*.nii'])
% delete([path 'forFSL/*.nii.gz'])
end


