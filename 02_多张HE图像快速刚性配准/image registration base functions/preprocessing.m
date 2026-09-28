function [im,impg,TA,fillval]=preprocessing(im,TA,szz,padall)
% pre-processing of H&E images for rigid registration.
% Written in 2020 by Ashley Lynn Kiemen, Johns Hopkins University
% please cite Kiemen et al, Nature methods (2022)
% last updated in December 2023 by ALK

    % 填色用背景像素中位数，避免 nested mode 把某通道取成 0
    fillval=background_fillval(im,TA);
    if ~isempty(padall)
        im=pad_im_both2(im,szz,padall,fillval);
        if size(TA,1)~=size(im,1) || size(TA,2)~=size(im,2)
            TA=pad_im_both2(TA,szz,padall,0);
        end
    end

    TA=TA>0;
    if ~isa(im,'uint8')
        im=im2uint8(im);
    end

    if size(im,3)==3
        ima=im(:,:,1);ima(~TA)=255;
        imb=im(:,:,2);imb(~TA)=255;
        imc=im(:,:,3);imc(~TA)=255;
        imp=cat(3,ima,imb,imc);
        impg=imcomplement(rgb2gray(imp));
    else
        imp=im;
        imp(~TA)=255;
        impg=imcomplement(imp);
    end

    impg=imgaussfilt(impg,2);
end


function fillval=background_fillval(im,TA)
    if size(im,3)==1
        im=repmat(im,[1 1 3]);
    end
    use_mask=exist('TA','var') && ~isempty(TA) && isequal(size(TA,1),size(im,1)) ...
        && isequal(size(TA,2),size(im,2)) && any(TA(:)==0);
    if use_mask
        bg=~(TA>0);
    else
        [h,w,~]=size(im);
        bw=max(8,round(min(h,w)*0.04));
        bg=false(h,w);
        bg(1:bw,:)=true; bg(end-bw+1:end,:)=true;
        bg(:,1:bw)=true; bg(:,end-bw+1:end)=true;
    end
    fillval=zeros(1,3);
    for c=1:3
        ch=im(:,:,c);
        fillval(c)=median(double(ch(bg)));
    end
end
