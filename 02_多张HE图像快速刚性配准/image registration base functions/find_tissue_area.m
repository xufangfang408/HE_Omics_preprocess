function [TA,fillval]=find_tissue_area(im0,ta_opt)
% Robust H&E tissue mask for rigid registration (coarse silhouette).
% 染色门限逐图 Otsu；保底和形态学参数由 A2 传入的 ta_opt 控制。
% Written in 2020 by Ashley Lynn Kiemen, Johns Hopkins University
% please cite Kiemen et al, Nature methods (2022)

    if ~exist('ta_opt','var') || isempty(ta_opt); ta_opt=struct(); end
    ta_opt=fill_ta_opt(ta_opt);

    if size(im0,3)==1
        im0=repmat(im0,[1 1 3]);
    end
    if ~isa(im0,'uint8')
        im0=im2uint8(im0);
    end
    [H,W,~]=size(im0);

    short_side=ta_opt.short_side;
    scale=short_side/min(H,W);
    if scale<1
        nh=max(32,round(H*scale));
        nw=max(32,round(W*scale));
        ims=im2double(imresize(im0,[nh nw],'bilinear'));
    else
        ims=im2double(im0);
        nh=H; nw=W;
    end
    ims=imgaussfilt(ims,ta_opt.gauss);

    bw=max(8,round(min(nh,nw)*0.06));
    border=cat(1, ...
        reshape(ims(1:bw,:,:),[],3), ...
        reshape(ims(end-bw+1:end,:,:),[],3), ...
        reshape(ims(:,1:bw,:),[],3), ...
        reshape(ims(:,end-bw+1:end,:),[],3));
    bg=median(border,1);
    fillval=uint8(round(bg*255));

    chroma=max(ims,[],3)-min(ims,[],3);
    dist=sqrt(sum((ims-reshape(double(bg),1,1,3)).^2,3));
    d01=mat2gray(dist);

    t_c=max(graythresh(chroma),ta_opt.chroma_min);
    t_d=max(graythresh(d01),ta_opt.dist_min);
    TA=(chroma>t_c) | (d01>t_d);

    r_close=max(3,round(min(nh,nw)/ta_opt.close_div));
    TA=imclose(TA,strel('disk',r_close));
    TA=imopen(TA,strel('disk',2));
    TA=fill_small_holes(TA,ta_opt.hole_frac);

    TA=bwareaopen(TA,max(16,round(0.002*nh*nw)));
    CC=bwconncomp(TA);
    if CC.NumObjects>0
        areas=cellfun(@numel,CC.PixelIdxList);
        keep=areas>=ta_opt.min_cc_frac*max(areas);
        TAnew=false(size(TA));
        TAnew(cat(1,CC.PixelIdxList{keep}))=true;
        TA=TAnew;
    end
    TA=fill_small_holes(TA,ta_opt.hole_frac);

    if scale<1
        TA=imresize(double(TA),[H W],'bilinear')>0.5;
        TA=imclose(TA,strel('disk',2));
        TA=fill_small_holes(TA,ta_opt.hole_frac);
    end

    fprintf('  tissue mask: %.1f%% of pixels (bg RGB = [%d %d %d])\n', ...
        100*mean(TA(:)), fillval(1), fillval(2), fillval(3));
    TA=uint8(TA);
end


function ta_opt=fill_ta_opt(ta_opt)
    if ~isfield(ta_opt,'chroma_min') || isempty(ta_opt.chroma_min); ta_opt.chroma_min=0.06; end
    if ~isfield(ta_opt,'dist_min') || isempty(ta_opt.dist_min); ta_opt.dist_min=0.12; end
    if ~isfield(ta_opt,'close_div') || isempty(ta_opt.close_div); ta_opt.close_div=80; end
    if ~isfield(ta_opt,'hole_frac') || isempty(ta_opt.hole_frac); ta_opt.hole_frac=0.02; end
    if ~isfield(ta_opt,'min_cc_frac') || isempty(ta_opt.min_cc_frac); ta_opt.min_cc_frac=0.04; end
    if ~isfield(ta_opt,'short_side') || isempty(ta_opt.short_side); ta_opt.short_side=768; end
    if ~isfield(ta_opt,'gauss') || isempty(ta_opt.gauss); ta_opt.gauss=1.5; end
end


function TA=fill_small_holes(TA, max_frac)
    filled=imfill(TA,'holes');
    holes=filled & ~TA;
    CC=bwconncomp(holes);
    if CC.NumObjects==0
        return;
    end
    tissue=max(sum(TA(:)),1);
    lim=max(32, max_frac*tissue);
    areas=cellfun(@numel,CC.PixelIdxList);
    for i=1:CC.NumObjects
        if areas(i)<=lim
            TA(CC.PixelIdxList{i})=true;
        end
    end
end
