function [im,TA]=get_ims(pth,nm,tp,recompute_ta,ta_opt,mask_file)
% loads RGB H&E image and loads or creates tissue mask.
% 优先使用 A2 降采样得到的 mask_file / pth/TA/<nm>.tif，与 HE 像素对齐。
% Written in 2020 by Ashley Lynn Kiemen, Johns Hopkins University
% please cite Kiemen et al, Nature methods (2022)
% last updated in December 2023 by ALK

    if ~exist('recompute_ta','var') || isempty(recompute_ta)
        recompute_ta=false;
    end
    if ~exist('ta_opt','var') || isempty(ta_opt)
        ta_opt=struct();
    end
    if ~exist([pth,'TA/'],'dir');mkdir([pth,'TA/']);end

    im=imread([pth,nm,tp]);
    if size(im,3)==1;im=cat(3,im,im,im);end
    if size(im,3)>3;im=im(:,:,1:3);end
    pthTA=[pth,'TA/'];
    tafile=[pthTA,nm,'.tif'];

    TA=[];
    if exist('mask_file','var') && ~isempty(mask_file) && exist(mask_file,'file') && ~recompute_ta
        TA=imread(mask_file);
    elseif (~recompute_ta) && exist(tafile,'file')
        TA=imread(tafile);
    end

    if isempty(TA)
        TA=find_tissue_area(im,ta_opt);
        imwrite(uint8(TA>0)*255, tafile);
    end
    if ndims(TA)==3
        TA=TA(:,:,1);
    end
    TA=uint8(TA>0);
    if ~isequal(size(TA), size(im(:,:,1)))
        error('mask 与 HE 尺寸不一致: %s%s  HE=[%d %d] mask=[%d %d]', ...
            nm, tp, size(im,1), size(im,2), size(TA,1), size(TA,2));
    end
end
