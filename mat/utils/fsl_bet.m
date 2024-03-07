function [ image_bet, mask_bet ] = fsl_bet(path,image,f,res)
%fsl bet2 function
%   filename : string filename WITHOUT nii, f: f value for bet
switch nargin
    case 1
        f = 0.5;
        res = [1,1,2];
    case 2
        res = [1,1,2];
end
[~, ~, ~] = mkdir([path 'forFSL']);
save_nii(make_nii(abs(rot90(image,-1)),res),[path 'forFSL/temp.nii'])
command = ['bet2 ' path 'forFSL/temp ' path 'forFSL/temp_bet_f' num2str(f) ' -m -f ' num2str(f)];
[status, cmdout] = system(command);
gunzip([path 'forFSL/temp_bet_f' num2str(f) '.nii.gz']);
gunzip([path 'forFSL/temp_bet_f' num2str(f) '_mask.nii.gz']);
imstr = load_nii([path 'forFSL/temp_bet_f' num2str(f) '.nii']);
image_bet = rot90(imstr.img);
maskstr = load_nii([path 'forFSL/temp_bet_f' num2str(f) '_mask.nii']);
mask_bet = double(rot90(maskstr.img));
delete([path 'forFSL/*.nii'])
delete([path 'forFSL/*.nii.gz'])
end
