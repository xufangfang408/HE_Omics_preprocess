function [R,rs,xy,amv]=group_of_reg(amv0,arf,iternum0,sz,rf,bb,opt)
% calculates sets of global registrations considering different initialization angles for a pair of greyscale images
% Written in 2020 by Ashley Lynn Kiemen, Johns Hopkins University
% please cite Kiemen et al, Nature methods (2022)
% last updated in December 2023 by ALK
%
% 评分用灰度相关 corr2（由 reg_ims_com 返回），不再覆盖成掩膜 IoU。
% 角度集合、提前停止阈值、分数地板由 opt 传入；缺省时保持可独立调用。

    if ~exist('bb','var') || isempty(bb);bb=0.9;end
    if ~exist('opt','var') || isempty(opt);opt=struct();end
    if ~isfield(opt,'init_angles') || isempty(opt.init_angles)
        opt.init_angles=[-2 177 87 268 -1 88 269 178 -7:2:7 179:183 89:93 270:272];
    end
    if ~isfield(opt,'score_floor') || isempty(opt.score_floor)
        opt.score_floor=-Inf;
    end
    if ~isfield(opt,'refine_theta') || isempty(opt.refine_theta)
        opt.refine_theta=-60:0.5:60;
    end
    if ~isfield(opt,'refine_thetaout') || isempty(opt.refine_thetaout)
        opt.refine_thetaout=2;
    end

    T=opt.init_angles;
    th={opt.refine_theta, opt.refine_thetaout};

    R=opt.score_floor;
    rs=0;
    xy=0;
    aa=arf==0;ab=amv0==0;
    arf=double(arf);amv0=double(amv0);
    arf=(arf-mean(arf(:)))/std(amv0(:));
    amv0=(amv0-mean(amv0(:)))/std(amv0(:));
    arf=arf-min(arf(:));arf(aa==1)=0;
    amv0=amv0-min(amv0(:));amv0(ab==1)=0;
    amv=amv0;
    for kp=1:length(T)
        try
            [~,amv1,rs1,xy1,RR]=reg_ims_com(amv0,arf,iternum0,sz,rf,T(kp),[0; 0],1,th);
            if isempty(RR) || (isscalar(RR) && isnan(RR));RR=opt.score_floor;end
        catch
            disp('catch')
            RR=opt.score_floor;
            amv1=amv0;rs1=0;xy1=[0; 0];
        end

        if RR>R;R=RR;rs=rs1;xy=xy1;amv=amv1;end
        if R>bb && kp>16;break;end
    end
end
