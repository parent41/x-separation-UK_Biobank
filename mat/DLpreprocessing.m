
function DLpreprocessing(path)

    % path = './tmp_chisep_2/';
    % path = './chisep_sub-1000011_ses-2/';

    addpath(genpath("./x-separation-UK_Biobank/"))
    addpath(genpath("./x-separation-UK_Biobank/mat/utils"))

    mkdir([path 'DLpreprocessing']);

    disp('Load data')

    voxelsize = [1.05 1 2];

    % Flip, register, and resample in bash

    % voxelsize_T1 = niftiinfo([path 't1_flair_mask/T1.nii.gz']).PixelDimensions;
    % voxelsize_Flair = niftiinfo([path 't1_flair_mask/T2_FLAIR.nii.gz']).PixelDimensions;
    % voxelsize_Mask = niftiinfo([path 't1_flair_mask/mask.nii.gz']).PixelDimensions;

    % T1_pre = fliplr( double(niftiread([path 't1_flair_mask/T1.nii.gz'])));
    % T1_pre = fliplr( double(niftiread([path 't1_flair_mask/T1.nii.gz'])));
    % Flair_pre = fliplr(double(niftiread([path 't1_flair_mask/T2_FLAIR.nii.gz'])));
    % Mask = fliplr( double(niftiread([path 't1_flair_mask/mask.nii.gz'])));

    % T1_pre = double(niftiread([path 't1_flair_mask/T1_flip.nii.gz']));
    % Flair_pre = double(niftiread([path 't1_flair_mask/T2_FLAIR_to_t1_flip.nii.gz']));
    % Mask = double(niftiread([path 't1_flair_mask/mask_flip.nii.gz']));

    % T1_pre_vox = resample_image(T1_pre,voxelsize_T1,voxelsize);
    % Flair_pre_vox = resample_image(Flair_pre,voxelsize_Flair,voxelsize);
    % Mask_pre_vox = resample_image(Mask,voxelsize_Mask,voxelsize);

    %%%% NEW: flip and resample outside of matlab

    T1_pre_vox = double(niftiread([path 't1_flair_mask/T1_flip_rsl.nii.gz']));
    Flair_pre_vox = double(niftiread([path 't1_flair_mask/T2_FLAIR_to_t1_flip_rsl.nii.gz']));
    Mask_pre_vox = double(niftiread([path 't1_flair_mask/mask_flip_rsl.nii.gz']));
    mask = Mask_pre_vox;

    % Added: extract brain
    T1w_bet_reg = T1_pre_vox.*Mask_pre_vox;
    Flair_bet_reg = Flair_pre_vox.*Mask_pre_vox;

    % [T1w_reg_bet, ~] = fsl_bet(path,T1_pre_vox(:,:,20:end),0.35,voxelsize);
    % [Flair_reg_bet, ~] = fsl_bet(path,Flair_pre_vox(:,:,20:end),0.35,voxelsize);
    % Flair_bet_reg = Flair_reg_bet;
    % [ T1w_bet_reg, ~] = fsl_flirt2(path,Flair_reg_bet,T1w_reg_bet,6,voxelsize);
    % matrix_size = [198,256,size(Flair_bet_reg,3)];
    % padsize = [198-size(T1w_bet_reg,1),256-size(T1w_bet_reg,2),0];

    % T1w_bet_reg = padarray(T1w_bet_reg,floor(padsize/2),0,'post');
    % Flair_bet_reg = padarray(Flair_bet_reg,floor(padsize/2),0,'post');
    % T1w_bet_reg = padarray(T1w_bet_reg,padsize-floor(padsize/2),0,'pre');
    % Flair_bet_reg = padarray(Flair_bet_reg,padsize-floor(padsize/2),0,'pre');

    mask = (T1w_bet_reg~=0).*(Flair_bet_reg~=0);

    nii_info = niftiinfo([path 'template_112_autocrop_new.nii.gz']);
    nii_info.Datatype = 'double';
    nii_info.ImageSize = size(T1w_bet_reg);
    nii_info.PixelDimensions = [1.05 1 2];
    niftiwrite(T1w_bet_reg, [path 'DLpreprocessing/T1w_bet_reg.nii'], nii_info, 'Compressed', true);
    niftiwrite(Flair_bet_reg, [path 'DLpreprocessing/Flair_bet_reg.nii'], nii_info, 'Compressed', true);

    % figure;imshow_3df(T1w_bet_reg.*mask,Flair_bet_reg.*mask)

    disp('Save data for DL of R2')

    mkdir([path 'DLpreprocessing'])
    save([path 'DLpreprocessing/registered_images.mat'],'T1w_bet_reg', 'Flair_bet_reg','mask');

    Mask = mask;
    t1_norm = T1w_bet_reg/prctile(T1w_bet_reg(T1w_bet_reg~=0),99).*Mask;
    flair_norm = Flair_bet_reg/prctile(Flair_bet_reg(Flair_bet_reg~=0),99).*Mask;

    r2_norm = Mask;
    % figure;imshow_3df(t1_norm,flair_norm,r2_norm,'range',[0 1]);

    niftiwrite(t1_norm, [path '/DLpreprocessing/t1_norm.nii'], nii_info, 'Compressed', true);
    niftiwrite(flair_norm, [path '/DLpreprocessing/flair_norm.nii'], nii_info, 'Compressed', true);

    save([path 'DLpreprocessing/Data_to_gpu_t2map.mat'],'r2_norm','flair_norm','t1_norm','Mask')
    
    %%% For T2star
    disp('Save data for DL of R2star')

    voxel_size = [0.8 0.8 3];
    voxelsize_new = [1.05 1 3];
    load([path 'GREprocessing/r2starimg_svd_2e.mat'])

    r2star_img = pd.*R2star_img;

    system(['bet2 ' path '/GREprocessing/pd.nii.gz ' path '/GREprocessing/pd -m -f 0.35']);
    mask = double(niftiread([path 'GREprocessing/pd_mask.nii.gz']));

    % save_nii(make_nii(abs(rot90(mean(R2star_img,4),-1)),voxelsize_new),[path 'DLpreprocessing/temp.nii'])
    % system(['bet2 ' path 'DLpreprocessing/temp.nii ' path 'DLpreprocessing/temp_bet -m -f 0.3']);
    % gunzip([path 'DLpreprocessing/temp_bet_mask.nii.gz']);
    % maskstr = load_nii([path 'DLpreprocessing/temp_bet_mask.nii']);
    % mask = double(rot90(maskstr.img));

    % [~,mask] = fsl_bet(mean(R2star_img,4),0.3,voxelsize_new);
    r2star_norm = mask;
    img = r2star_img.*mask/prctile(r2star_img(mask~=0),99);

    Mask = mask;
    save([path 'DLpreprocessing/Data_to_gpu_t2starmap.mat'],'r2star_norm','img','Mask')
    % figure;imshow_3df(img)
    % mkdir([path 'DLpreprocessing/inf_from_gpu'])

    nii_info = niftiinfo([path 'template_113_autocrop_new.nii.gz']);
    nii_info.Datatype = 'double';
    nii_info.ImageSize = size(r2star_norm);
    nii_info.PixelDimensions = [1.05 1 3];
    niftiwrite(r2star_norm, [path '/DLpreprocessing/r2star_norm.nii'], nii_info, 'Compressed', true);
    niftiwrite(Mask, [path '/DLpreprocessing/Mask.nii'], nii_info, 'Compressed', true);

end

