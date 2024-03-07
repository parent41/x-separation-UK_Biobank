
function Chisep(path,xfm)

    % path = './tmp_chisep_2/';
    % path = 'chisep_sub-1000011_ses-2/';
    % xfm = './preproc/swi/sub-1000011_ses-2_SWI_linregT1_0_GenericAffine.xfm';

    addpath(genpath("./x-separation-UK_Biobank/"))
    addpath(genpath("./x-separation-UK_Biobank/mat/utils"))

    disp('Load data')

    load([path 'GREprocessing/for_xsep.mat'],'dB_vsf','r2star','Mag_brain','mask_vsf','TE','B0','H','voxsz','voxelsize_new','voxel_size')
    load([path 'GREprocessing/QSM.mat'],'x_sa')
    params.voxel_size = voxsz; % [mm unit]
    params.CF = 123e6; 
    params.b0_dir = H; 
    params.Dr_pos = 137; 
    params.Dr_neg = params.Dr_pos;
    params.lambda = 1; 
    params.lambda_CSF = 1; 

    mkdir([path 'reconstruction'])

    %% T2: MPRAGE+FLAIR Ds, T2*: DL Ds

    disp('Combine data')

    dire = dir([path 'r2_r2star_mapping/inf*']);
    load([path 'r2_r2star_mapping/' dire(1).name])
    load([path 'r2_r2star_mapping/' dire(2).name])
    r2star_deep = permute(-10*(log(R2star+eps)),[2 3 1]).*mask_vsf;
    load([path 'DLpreprocessing/registered_images.mat'])
    R2_MPR_FLA = -10*log(permute(T2,[2 3 1]));
    R2_MPR_FLA(isinf(R2_MPR_FLA) | isnan(R2_MPR_FLA))=0;
    % R2_MPR_FLA_resample = resample_image(rot90(R2_MPR_FLA,-2),voxel_size,voxelsize_new);
    % R2_MPR_FLA_resample = R2_MPR_FLA;

    nii_info = niftiinfo([path 'template_112_autocrop_new.nii.gz']);
    nii_info.Datatype = 'double';
    nii_info.ImageSize = size(R2_MPR_FLA);
    nii_info.PixelDimensions = [1.05 1 2];
    niftiwrite(R2_MPR_FLA, [path '/reconstruction/R2_MPR_FLA.nii'], nii_info, 'Compressed', true);

    system(['nii2mnc -clobber ' path '/reconstruction/R2_MPR_FLA.nii.gz ' path '/reconstruction/R2_MPR_FLA.mnc']);
    system(['autocrop -clobber -double -step 1.05 1 3 ' path '/reconstruction/R2_MPR_FLA.mnc ' path '/reconstruction/R2_MPR_FLA_resample.mnc']);
    system(['mnc2nii ' path '/reconstruction/R2_MPR_FLA_resample.mnc ' path '/reconstruction/R2_MPR_FLA_resample.nii']);
    system(['gzip -f ' path '/reconstruction/R2_MPR_FLA_resample.nii']);

    R2_MPR_FLA_resample = niftiread([path '/reconstruction/R2_MPR_FLA_resample.nii.gz']);

    nii_info = niftiinfo([path 'template_113_autocrop_new.nii.gz']);
    nii_info.Datatype = 'double';
    nii_info.ImageSize = size(r2star_deep);
    nii_info.PixelDimensions = [1.05 1 3];
    niftiwrite(r2star_deep, [path '/reconstruction/r2star_deep.nii'], nii_info, 'Compressed', true);

    % system(['antsApplyTransforms -i ' path '/reconstruction/R2_MPR_FLA.nii.gz -r ' path '/reconstruction/r2star_deep.nii.gz -o ' path '/reconstruction/r2_MPR_FLA_reg_r2s.nii.gz -t [ ' xfm ',1 ]']);

    % nii_info.ImageSize = size(R2_MPR_FLA_resample);
    % niftiwrite(R2_MPR_FLA_resample, [path '/reconstruction/R2_MPR_FLA_resample.nii'], nii_info, 'Compressed', true);

    % [ r2_MPR_FLA_reg_r2s, ~] = fsl_flirt2(path,r2star_deep.*mask_vsf,R2_MPR_FLA_resample,6,voxelsize_new);
    % r2_MPR_FLA_reg_r2s = double(r2_MPR_FLA_reg_r2s);
    % [ r2_MPR_FLA_reg_r2s, ~] = fsl_flirt2( r2star_deep.*mask_vsf,rot90(R2_MPR_FLA,-2),6,voxelsize_new);
    % system(['flirt -in ' path '/reconstruction/R2_MPR_FLA_resample.nii.gz -ref ' path '/reconstruction/r2star_deep.nii.gz -dof 6 -out ' path '/reconstruction/r2_MPR_FLA_reg_r2s_new.nii.gz  -interp sinc']);
    system(['antsApplyTransforms -i ' path '/reconstruction/R2_MPR_FLA_resample.nii.gz -r ' path '/reconstruction/r2star_deep.nii.gz -o ' path '/reconstruction/r2_MPR_FLA_reg_r2s_new.nii.gz -t [ ' path '/swi_to_t1_ANTs_t1flip_0GenericAffine.mat,1 ]']);
    r2_MPR_FLA_reg_r2s = niftiread([path '/reconstruction/r2_MPR_FLA_reg_r2s_new.nii.gz']);
    r2prime_MPR_FLA_r2sdeep = (r2star_deep-r2_MPR_FLA_reg_r2s).*mask_vsf;
    r2prime_MPR_FLA_r2sdeep(r2prime_MPR_FLA_r2sdeep<0) = 0;

    nii_info = niftiinfo([path 'template_113_autocrop_new.nii.gz']);
    nii_info.Datatype = 'double';
    nii_info.ImageSize = size(r2_MPR_FLA_reg_r2s);
    nii_info.PixelDimensions = [1.05 1 3];
    % niftiwrite(r2_MPR_FLA_reg_r2s, [path '/reconstruction/r2_MPR_FLA_reg_r2s.nii'], nii_info, 'Compressed', true);
    niftiwrite(r2prime_MPR_FLA_r2sdeep, [path '/reconstruction/r2prime_MPR_FLA_r2sdeep.nii'], nii_info, 'Compressed', true);

    disp('Chi separation')

    [x_pos_MPR_FLA_r2sdeep, x_neg_MPR_FLA_r2sdeep, x_tot_MPR_FLA_r2sdeep]= x_sep_SA(dB_vsf/2/pi/TE,r2prime_MPR_FLA_r2sdeep, mask_vsf,params,x_sa);
    % figure(1);imshow_3df(fliplr(x_pos_MPR_FLA_r2sdeep),[0 .1], fliplr(x_neg_MPR_FLA_r2sdeep),[0 .1],fliplr(x_tot_MPR_FLA_r2sdeep),[-.1 .1])
    save([path 'reconstruction/x_sep_MPRAGE_FLAIR_r2stardeep.mat'],'r2star_deep','r2_MPR_FLA_reg_r2s','r2prime_MPR_FLA_r2sdeep','x_pos_MPR_FLA_r2sdeep','x_neg_MPR_FLA_r2sdeep','x_tot_MPR_FLA_r2sdeep')
    clear R2 R2_MPR_FLA r2_MPR_FLA_reg_r2s R2_MPR_FLA_resample r2_zpd r2prime_MPR_FLA_r2sdeep R_mat2 T1w_bet_reg T2 t2w_reg t2w_resample t2w_zpd  x_neg_ds_r2s x_pos_ds_r2s x_tot_ds_r2s mask Flair_bet_reg R2star r2star_deep

    nii_info.ImageSize = size(x_pos_MPR_FLA_r2sdeep);
    niftiwrite(x_pos_MPR_FLA_r2sdeep, [path '/reconstruction/x_pos_MPR_FLA_r2sdeep.nii'], nii_info, 'Compressed', true);
    niftiwrite(x_neg_MPR_FLA_r2sdeep, [path '/reconstruction/x_neg_MPR_FLA_r2sdeep.nii'], nii_info, 'Compressed', true);
    niftiwrite(x_tot_MPR_FLA_r2sdeep, [path '/reconstruction/x_tot_MPR_FLA_r2sdeep.nii'], nii_info, 'Compressed', true);


    % % Save to niftii
    % load([path 'reconstruction/x_sep_MPRAGE_FLAIR_r2stardeep.mat'],'r2star_deep','r2_MPR_FLA_reg_r2s','r2prime_MPR_FLA_r2sdeep','x_pos_MPR_FLA_r2sdeep','x_neg_MPR_FLA_r2sdeep','x_tot_MPR_FLA_r2sdeep')
    % save_nii(make_nii(r2star_deep,voxelsize_new),[path 'reconstruction/r2star_deep.nii'])
    % save_nii(make_nii(r2_MPR_FLA_reg_r2s,voxelsize_new),[path 'reconstruction/r2_MPR_FLA_reg_r2s.nii'])
    % save_nii(make_nii(r2prime_MPR_FLA_r2sdeep,voxelsize_new),[path 'reconstruction/r2prime_MPR_FLA_r2sdeep.nii'])
    % save_nii(make_nii(x_pos_MPR_FLA_r2sdeep,voxelsize_new),[path 'reconstruction/x_pos_MPR_FLA_r2sdeep.nii'])
    % save_nii(make_nii(x_neg_MPR_FLA_r2sdeep,voxelsize_new),[path 'reconstruction/x_neg_MPR_FLA_r2sdeep.nii'])
    % save_nii(make_nii(x_tot_MPR_FLA_r2sdeep,voxelsize_new),[path 'reconstruction/x_tot_MPR_FLA_r2sdeep.nii'])


    %% T2: Cons, T2*: DL

    % load([path 'r2_r2star_mapping/' dire(2).name])
    % r2star_deep = permute(-10*(log(R2star+eps)),[2 3 1]).*mask_vsf;

    % r2_con = 13.98;
    % r2prime_con = r2star_deep.*mask_vsf-r2_con;
    % r2prime_con(r2prime_con<0) = 0;
    % [x_pos_sa_con, x_neg_sa_con, x_tot_sa_con]= x_sep_SA(dB_vsf/2/pi/TE,r2prime_con, mask_vsf,params,x_sa);
    % % figure(2);imshow_3df(fliplr(x_pos_sa_con),[0 .1], fliplr(x_neg_sa_con),[0 .1],fliplr(x_tot_sa_con),[-.1 .1])

    % save([path 'reconstruction/x_sep_constant_r2stardeep_cal.mat'],'r2_con','r2prime_con','x_pos_sa_con','x_neg_sa_con','x_tot_sa_con')

end


