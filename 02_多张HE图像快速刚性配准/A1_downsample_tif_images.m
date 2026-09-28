function A1_downsample_tif_images(he_dir, mask_dir, downsample_name, downsample_factor, redo)
% 已并入 A2_calculate_image_registration。
% 单独跑降采样时仍可用本函数；完整流程请直接运行 A2（会先降采样再配准）。
%
% HE 与 mask 使用同一目标尺寸、最近邻，保证降采样后组织区域对齐。

    if nargin < 1 || isempty(he_dir)
        he_dir = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/interim/PGD_dataset/HE_unregistered";
    end
    if nargin < 2 || isempty(mask_dir)
        mask_dir = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/interim/PGD_dataset/HE_unregistered_masks";
    end
    if nargin < 3 || isempty(downsample_name)
        downsample_name = "downsample_image";
    end
    if nargin < 4 || isempty(downsample_factor)
        downsample_factor = 2;
    end
    if nargin < 5 || isempty(redo)
        redo = false;
    end

    path(path,'/mnt/zzf_nas/A_xff/Spaceland-omics/spaceland_omics/rigid_registration/image registration base functions');
    downsample_he_and_masks(he_dir, mask_dir, downsample_name, downsample_factor, redo);
end
