function [imout,tform,cent,f,Rout,info]=calculate_global_reg(imrf,immv,rf,iternum,opt)
% calculates global registration of a pair of greyscale, downsampled H&E images.
% imrf == reference image
% immv == moving image
% rf == reduce images by _ times
% iternum == number of iterations of registration code
% opt == 由 A2 传入的参数结构体；缺省字段在本函数内补齐
% Written in 2020 by Ashley Lynn Kiemen, Johns Hopkins University
% please cite Kiemen et al, Nature methods (2022)

    if ~exist('opt','var') || isempty(opt);opt=struct();end
    opt=fill_missing_reg_opt(opt);
    bb=opt.early_stop_bb;

    % pre-registration image processing
    amv=imresize(immv,1/rf);amv=imgaussfilt(amv,2);
    arf=imresize(imrf,1/rf);arf=imgaussfilt(arf,2);
    sz=[0 0];cent=[0 0];

    iternum0=opt.iternum0;
    nfine=max(iternum-iternum0,1);
    th={opt.refine_theta, opt.refine_thetaout};

    % 未翻转：粗配准 + 精配准
    [Ccoarse1,rs1,xy1]=group_of_reg(amv,arf,iternum0,sz,rf,bb,opt);
    [tform1,amvout1,~,~,Cfine1]=reg_ims_com(amv,arf,nfine,sz,rf,rs1,xy1,0,th);
    iou1=mask_iou(arf,amvout1);

    Ccoarse2=NaN;Cfine2=NaN;iou2=NaN;
    f=0;
    tform=tform1;amvout=amvout1;rs=rs1;xy=xy1;
    Rout=Cfine1;

    if opt.always_try_flip
        amv2=amv(end:-1:1,:,:);
        [Ccoarse2,rs2,xy2]=group_of_reg(amv2,arf,iternum0,sz,rf,bb,opt);
        [tform2,amvout2,~,~,Cfine2]=reg_ims_com(amv2,arf,nfine,sz,rf,rs2,xy2,0,th);
        iou2=mask_iou(arf,amvout2);

        % 用精配准灰度相关选翻转，掩膜 IoU 只打印不决策
        if Cfine2 > Cfine1 + opt.flip_score_margin
            f=1;
            tform=tform2;amvout=amvout2;rs=rs2;xy=xy2;amv=amv2;
            Rout=Cfine2;
        end
    end

    fprintf(['  flip compare  coarse corr: unflip=%.3f  flip=%.3f\n' ...
             '               fine   corr: unflip=%.3f  flip=%.3f  margin=%.3f  -> f=%d\n' ...
             '               fine   IoU : unflip=%.3f  flip=%.3f\n'], ...
        Ccoarse1,Ccoarse2,Cfine1,Cfine2,opt.flip_score_margin,f,iou1,iou2);

    % create output image
    Rin=imref2d(size(immv));
    if sum(abs(cent))==0
      mx=mean(Rin.XWorldLimits);
      my=mean(Rin.YWorldLimits);
      cent=[mx my];
    end
    Rin.XWorldLimits = Rin.XWorldLimits-cent(1);
    Rin.YWorldLimits = Rin.YWorldLimits-cent(2);

    if f==1
        immv=immv(end:-1:1,:,:);
    end

    imout=imwarp(immv,Rin,tform,'nearest','Outputview',Rin,'Fillvalues',0);

    info=struct();
    info.f=f;
    info.corr_coarse_unflip=Ccoarse1;
    info.corr_coarse_flip=Ccoarse2;
    info.corr_fine_unflip=Cfine1;
    info.corr_fine_flip=Cfine2;
    info.corr_chosen=Rout;
    info.iou_fine_unflip=iou1;
    info.iou_fine_flip=iou2;
    info.rs=rs;
    info.xy=xy;
end


function opt=fill_missing_reg_opt(opt)
    if ~isfield(opt,'iternum0') || isempty(opt.iternum0);opt.iternum0=2;end
    if ~isfield(opt,'always_try_flip') || isempty(opt.always_try_flip);opt.always_try_flip=true;end
    if ~isfield(opt,'flip_score_margin') || isempty(opt.flip_score_margin);opt.flip_score_margin=0.02;end
    if ~isfield(opt,'early_stop_bb') || isempty(opt.early_stop_bb);opt.early_stop_bb=0.9;end
    if ~isfield(opt,'refine_theta') || isempty(opt.refine_theta);opt.refine_theta=-60:0.5:60;end
    if ~isfield(opt,'refine_thetaout') || isempty(opt.refine_thetaout);opt.refine_thetaout=2;end
    if ~isfield(opt,'init_angles') || isempty(opt.init_angles)
        opt.init_angles=[-2 177 87 268 -1 88 269 178 -7:2:7 179:183 89:93 270:272];
    end
    if ~isfield(opt,'score_floor') || isempty(opt.score_floor);opt.score_floor=-Inf;end
end


function iou=mask_iou(a,b)
    aa=double(a>0)+double(b>0);
    un=sum(aa(:)>0);
    if un==0;iou=0;else;iou=sum(aa(:)==2)/un;end
end
