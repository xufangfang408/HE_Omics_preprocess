function A2_calculate_image_registration(he_dir_in, mask_dir_in, ds_name_in, ds_factor_in, zc_in)

% 2025/10/28 更新：添加叠加可视化模块和注册评估指标计算模块
% 2025/12/10 更新：添加IoU对比模块，如果IoU_after < IoU_before，则不使用配准结果，保留原始图像
% 2025/12/10 更新：删掉所有和弹性配准相关的代码
% 2025/12/16 更新：保存的overlay叠加图像改为uint8格式，减少存储空间
% 2026/3/10 更新：只保留IoU评估指标，删除其余指标
% 2026/3/10 更新：删除输出数据结构，去掉所有和弹性配准相关的代码
% 2026/9/2  更新：翻转始终评估；用灰度相关选翻转；参数全部集中到本文件开头；修复 overlay 背景发绿
% 2026/9/26 更新：A1 降采样并入本文件；HE 与 mask 用同一尺寸最近邻降采样，配准时同步 pad/flip/warp

%% =====================================================================
%  一、必填输入（数据从哪来、降采样到哪；与下面配准调参分开）
%  函数入参非空时覆盖对应项。
%  =====================================================================

    % HE_dir: 原分辨率 H&E 目录（只读顶层 *.tif）
    HE_dir = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/interim/PGD_dataset/HE_unregistered";
    % mask_dir: 与 HE 同名的组织掩膜目录。会依次找 mask_dir/<fname> 和 mask_dir/mask/<fname>
    mask_dir = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/interim/PGD_dataset/HE_unregistered_masks";
    % downsample_name: 降采样结果子文件夹名。HE 写到 HE_dir/<name>/，mask 写到 mask_dir/<name>/
    downsample_name = "downsample_image";
    % downsample_factor: 边长缩小倍数。2 = 长宽各 /2。HE 与 mask 共用同一目标尺寸
    downsample_factor = 2;
    % redo_downsample: true=覆盖已有降采样文件；false=已存在则跳过
    redo_downsample = false;
    % zc: 参考切片在文件列表中的序号；[] = 自动取中间一张
    zc = 47;

%% =====================================================================
%  二、配准调参（只改算法行为，不必动上面的路径）
%  =====================================================================

    % ----- 预处理 -----
    % padall: 每张图四周填充像素（降采样分辨率），避免旋转/平移后组织出界
    padall = 250;
    % recompute_ta: true=忽略已提供的 mask，改用 find_tissue_area 重算。
    %              默认 false：配准使用与 HE 同步降采样后的外部 mask。
    recompute_ta = false;

    % ----- 组织掩膜兜底（仅当某张图缺少外部 mask，或 recompute_ta=true） -----
    % 只在「整批系统性偏差」时改，改完保持 recompute_ta=true 重跑。不必按张手调。
    %
    % 浅染伊红大面积缺失（组织被当成背景）
    %     → 略降 chroma_min（如 0.06→0.04）、dist_min（如 0.12→0.08）
    % 灰底 / 扫描棋盘格被当成组织
    %     → 略升 chroma_min（如 0.06→0.08）、dist_min（如 0.12→0.16）
    % C 形凹槽被填实、掩膜变成粗多边形
    %     → 增大 close_div（如 80→120，闭运算更小）或降低 hole_frac（如 0.02→0.01）
    % 小组织块 / 游离碎片被丢掉
    %     → 降低 min_cc_frac（如 0.04→0.02）
    % 掩膜边缘呈方块、贴边不准
    %     → 略增 ta_short_side（如 768→1024，更慢）
    chroma_min   = 0.06;  % Otsu(chroma) 的下限。越低越容易把浅染算进组织
    dist_min     = 0.12;  % Otsu(相对背景色差) 的下限。越高越能挡住灰底噪声
    close_div    = 80;    % 闭运算半径 = min(高,宽)/close_div。越大半径越小，越不容易桥接 C 形开口
    hole_frac    = 0.02;  % 只填面积 < 组织面积×该比例 的封闭孔。过大则 C 形凹槽也会被填
    min_cc_frac  = 0.04;  % 丢掉面积 < 最大连通域×该比例 的碎块
    ta_short_side = 768;  % 掩膜计算短边像素。越大越贴边、越慢
    ta_gauss      = 1.5;  % 掩膜前高斯平滑。棋盘格重可略增（如 2～3）

    % ----- 多尺度配准 -----
    % rsc: 配准在 1/rsc 分辨率上进行。越大越快，但翻转/180° 判决越不准。
    %      原值 9；翻转经常判错时改为 6 或 4。
    rsc = 6;
    % iternum: 粗+精配准的总迭代次数
    iternum = 6;
    % iternum0: 粗配准迭代。过小（如 1）会让翻转分支的粗分数不可靠。建议 2。
    iternum0 = 2;

    % ----- 角度搜索 -----
    % init_angles: 粗配准的初始旋转（度）。必须包含 ~90/180/270，否则
    %              「上下翻转 + 180°」这条左右镜像路径走不到。
    init_angles = [-2 177 87 268 -1 88 269 178 -7:2:7 179:183 89:93 270:272];
    % refine_theta: 精配准相对当前角度的搜索范围，不是绝对角。
    %               大角度由 init_angles 提供，这里只做 ±60° 微调。
    refine_theta = -60:0.5:60;
    refine_thetaout = 2;

    % ----- 翻转判决 -----
    % always_try_flip: true=每对图都跑「不翻转」和「上下翻转」两套精配准，
    %                  再用灰度相关选 f。不要用掩膜 IoU 门槛决定「要不要试翻转」。
    always_try_flip = true;
    % flip_score_margin: 翻转的精配准 corr2 必须比不翻转高这么多才采用，
    %                    避免噪声导致误翻。真正镜像时分差通常远大于 0.02。
    flip_score_margin = 0.02;
    % score_floor: 粗配准角度搜索的初始分数。原代码用 0.2，真实 corr 低于
    %              0.2 时分数被地板盖住，R2>R 永远不成立。必须用 -Inf。
    score_floor = -Inf;
    % early_stop_bb: 某个初始角的 corr 超过此值且已试过足够角度则提前停。
    early_stop_bb = 0.9;

    % ----- 参考图 A/B/C 切换 -----
    % ref_switch_ct: 精配准灰度相关低于此值，则再试更早的参考图。
    %                现在比较的是 corr2，不再是掩膜 IoU；0.9 对 corr 过严，
    %                会频繁试 B/C。建议 0.6。
    ref_switch_ct = 0.6;

    % ----- 配准验收（掩膜 IoU，只用于「要不要丢掉这次结果」） -----
    % iou_min_abs: 配准后 IoU 绝对下限。只要求「比配准前好」不够，
    %              镜像未翻时轮廓对上也能从 0.24 升到 0.46。低于此值回退。
    iou_min_abs = 0.50;
    % iou_require_improve: true 时仍要求 IoU_after >= IoU_before
    iou_require_improve = true;

    % ----- 输出与可视化 -----
    tpout = '.jpg';
    % overlay 通道强度。必须 <=1，且配准前/后都先把背景置 0、组织亮度拉齐，
    % 否则白色画布会灌进绿色通道，整张 overlay 发绿。
    overlay_alpha_ref = 1.0;   % 参考图 → 红
    overlay_alpha_mov = 1.0;   % 移动图 → 绿；重叠 → 黄

%% =====================================================================
%  把上面的参数打进 opt，后面和子函数都只读 opt
%  =====================================================================
    if nargin>=1 && ~isempty(he_dir_in); HE_dir=he_dir_in; end
    if nargin>=2 && ~isempty(mask_dir_in); mask_dir=mask_dir_in; end
    if nargin>=3 && ~isempty(ds_name_in); downsample_name=ds_name_in; end
    if nargin>=4 && ~isempty(ds_factor_in); downsample_factor=ds_factor_in; end
    if nargin>=5 && ~isempty(zc_in); zc=zc_in; end

    opt = struct();
    opt.iternum0 = iternum0;
    opt.init_angles = init_angles;
    opt.refine_theta = refine_theta;
    opt.refine_thetaout = refine_thetaout;
    opt.always_try_flip = always_try_flip;
    opt.flip_score_margin = flip_score_margin;
    opt.score_floor = score_floor;
    opt.early_stop_bb = early_stop_bb;
    opt.ta = struct();
    opt.ta.chroma_min = chroma_min;
    opt.ta.dist_min = dist_min;
    opt.ta.close_div = close_div;
    opt.ta.hole_frac = hole_frac;
    opt.ta.min_cc_frac = min_cc_frac;
    opt.ta.short_side = ta_short_side;
    opt.ta.gauss = ta_gauss;

    path(path,'/mnt/zzf_nas/A_xff/Spaceland-omics/spaceland_omics/rigid_registration/image registration base functions');

    warning('off','all');
    [pth, mask_ds_dir] = downsample_he_and_masks( ...
        HE_dir, mask_dir, downsample_name, downsample_factor, redo_downsample);
    if pth(end)~='/';pth=[pth,'/'];end
    if mask_ds_dir(end)~='/';mask_ds_dir=[mask_ds_dir,'/'];end

    imlist=dir([pth,'*tif']);
    if isempty(imlist);imlist=dir([pth,'*jpg']);end
    if isempty(imlist);disp('no images found');return;end
    tp=imlist(1).name(end-3:end);

    if ~exist('zc','var') || isempty(zc);zc=ceil(length(imlist)/2);end
    rf=[zc:-1:2 zc:length(imlist)-1 0];
    mv=[zc-1:-1:1 zc+1:length(imlist)];

    szz=[0 0];
    for kk=1:length(imlist)
        inf=imfinfo([pth,imlist(kk).name]);
        szz=[max([szz(1),inf.Height]) max([szz(2),inf.Width])];
    end

    outpthG=[pth,'registered/'];
    matpth=[outpthG,'save_warps/'];
    outTA=[outpthG,'TA/'];
    mkdir(outpthG);mkdir(matpth);mkdir(outTA);

    nm=imlist(zc).name(1:end-4);
    [imzc,TAzc]=get_ims(pth,nm,tp,recompute_ta,opt.ta,[mask_ds_dir,nm,'.tif']);
    [imzc,imzcg,TAzc]=preprocessing(imzc,TAzc,szz,padall);
    disp(['Reference image: ',nm])

    imwrite(imzc,[outpthG,nm,tpout]);
    imwrite(uint8(TAzc>0)*255,[outTA,nm,'.tif']);
    save([matpth,nm,'.mat'],'zc','szz','padall','downsample_factor');

    img=imzcg;TA=TAzc;
    img0=imzcg;TA0=TAzc;krf0=zc;
    img00=imzcg;TA00=TAzc;krf00=zc;

    metrics = repmat(struct('mov','','ref','','f',NaN, ...
        'corr_unflip',NaN,'corr_flip',NaN,'corr_chosen',NaN, ...
        'IoU_before',NaN,'IoU_after',NaN,'reverted',false), length(mv), 1);

    for kk=1:length(mv)
        t1=tic;
        fprintf(['Image ',num2str(kk),' of ',num2str(length(imlist)-1),...
            '\n  reference image:  ',imlist(rf(kk)).name(1:end-4),...
            '\n  moving image:  ',imlist(mv(kk)).name(1:end-4),'\n']);
        nm=imlist(mv(kk)).name(1:end-4);
        [immv0,TAmv]=get_ims(pth,nm,tp,recompute_ta,opt.ta,[mask_ds_dir,nm,'.tif']);
        [immv,immvg,TAmv,fillval]=preprocessing(immv0,TAmv,szz,padall);

        if rf(kk)==zc
            imrfgA=img;TArfA=TA;krfA=zc;
            imrfgB=img0;TArfB=TA0;krfB=krf0;
            imrfgC=img00;TArfC=TA00;krfC=krf00;
        end

        if exist([matpth,'D/',nm,'.mat'],'file')
            disp('Registration already calculated');disp('please delete registered folder to recalculate')
        else
            RB=0.4;RC=0.4;immvGgB=immvg;immvGgC=immvg;
            infoB=struct('corr_fine_unflip',NaN,'corr_fine_flip',NaN,'corr_chosen',NaN,'f',0);
            infoC=infoB;

            [immvGg,tform,cent,f,R,info]=calculate_global_reg(imrfgA,immvg,rsc,iternum,opt);
            if R<ref_switch_ct
                [immvGgB,tformB,centB,fB,RB,infoB]=calculate_global_reg(imrfgB,immvg,rsc,iternum,opt);
                disp('RB');
            end
            if R<ref_switch_ct && RB<ref_switch_ct
                [immvGgC,tformC,centC,fC,RC,infoC]=calculate_global_reg(imrfgC,immvg,rsc,iternum,opt);
                disp('RC');
            end

            RR=[R RB RC];
            [~,ii]=max(RR);disp(RR)
            if ii==1
                imrfg=imrfgA;TArf=TArfA;krf=krfA;disp('chose image A')
                offset=0;
            elseif ii==2
                immvGg=immvGgB;tform=tformB;cent=centB;f=fB;info=infoB;disp('chose image B')
                imrfg=imrfgB;TArf=TArfB;krf=krfB;
                offset=1;
            else
                immvGg=immvGgC;tform=tformC;cent=centC;f=fC;info=infoC;disp('chose image C')
                imrfg=imrfgC;TArf=TArfC;krf=krfC;
                offset=2;
            end

            fprintf('  chosen f=%d  corr_unflip=%.3f  corr_flip=%.3f  corr_chosen=%.3f\n', ...
                f, info.corr_fine_unflip, info.corr_fine_flip, info.corr_chosen);

            immvG=register_global_im(immv,tform,cent,f,fillval);
            TAmvG=register_global_im(TAmv,tform,cent,f,0);

            [IoU_before, IoU_after, ref_gray, mov_gray, reg_gray] = ...
                evaluate_registration(imrfg, immv, immvG, TArf, TAmv, TAmvG);

            reverted = false;
            fail_improve = iou_require_improve && (IoU_after < IoU_before);
            fail_abs = IoU_after < iou_min_abs;
            if fail_improve || fail_abs
                if fail_abs && ~fail_improve
                    disp('IoU_after below absolute floor, revert to unregistered image');
                elseif fail_improve && fail_abs
                    disp('IoU_after < IoU_before and below floor, revert to unregistered image');
                else
                    disp('IoU_after < IoU_before, revert to unregistered image');
                end
                immvG = immv;
                TAmvG = TAmv;
                immvGg = immvg;
                tform = affine2d(eye(3));
                cent = [0 0];
                f = 0;
                reverted = true;
                [IoU_before, IoU_after, ref_gray, mov_gray, reg_gray] = ...
                    evaluate_registration(imrfg, immv, immvG, TArf, TAmv, TAmvG);
            end

            save([matpth,nm,'.mat'],'tform','f','cent','szz','padall','krf','info','downsample_factor');
            imwrite(immvG,[outpthG,nm,tpout]);
            imwrite(uint8(TAmvG>0)*255,[outTA,nm,'.tif']);

            metrics(kk).mov = imlist(mv(kk)).name;
            metrics(kk).ref = imlist(rf(kk-offset)).name;
            metrics(kk).f = f;
            metrics(kk).corr_unflip = info.corr_fine_unflip;
            metrics(kk).corr_flip = info.corr_fine_flip;
            metrics(kk).corr_chosen = info.corr_chosen;
            metrics(kk).IoU_before = IoU_before;
            metrics(kk).IoU_after  = IoU_after;
            metrics(kk).reverted = reverted;

            fprintf('IoU before: %.3f, after: %.3f  (floor=%.2f, reverted=%d)\n', ...
                IoU_before, IoU_after, iou_min_abs, reverted);

            overlay_dir = fullfile(outpthG, 'overlay');
            if ~exist(overlay_dir, 'dir'); mkdir(overlay_dir); end

            overlay_before = compose_overlay(ref_gray, mov_gray, overlay_alpha_ref, overlay_alpha_mov);
            overlay_after  = compose_overlay(ref_gray, reg_gray, overlay_alpha_ref, overlay_alpha_mov);
            combined = cat(2, overlay_before, overlay_after);
            combined_uint8 = im2uint8(min(max(combined,0),1));

            out_name = sprintf('overlay_%s_vs_%s.jpg', imlist(mv(kk)).name(1:end-4), imlist(rf(kk-offset)).name(1:end-4));
            imwrite(combined_uint8, fullfile(overlay_dir, out_name));
            fprintf('Saved: %s\n', fullfile(overlay_dir, out_name));
        end

        imrfgC=imrfgB;TArfC=TArfB;krfC=krfB;
        imrfgB=imrfgA;TArfB=TArfA;krfB=krfA;
        imrfgA=immvGg;TArfA=TAmvG;krfA=mv(kk);
        if mv(kk)==mv(1);img0=immvGg;TA0=TAmvG;krf0=mv(kk);end
        try if mv(kk)==mv(2);img00=immvGg;TA00=TAmvG;krf00=mv(kk);end;catch;end

        out_table = struct2table(metrics);
        writetable(out_table, fullfile(outpthG, 'registration_metrics.csv'));
        disp(['Registration metrics saved to ', fullfile(outpthG, 'registration_metrics.csv')]);

        toc(t1);
    end
    warning('on','all');
end


function [IoU_before, IoU_after, ref_gray, mov_gray, reg_gray] = ...
        evaluate_registration(imrfg, immv, immvG, TArf, TAmv, TAmvG)

    ref_mask = TArf ~= 0;
    mov_mask = TAmv ~= 0;
    reg_mask = TAmvG ~= 0;

    % 配准前/后的灰度都按「组织亮、背景 0」统一，避免 RGB 白底灌进 overlay
    ref_gray = gray_for_overlay(imrfg, ref_mask);
    mov_gray = gray_for_overlay(immv,  mov_mask);
    reg_gray = gray_for_overlay(immvG, reg_mask);

    union_before = ref_mask | mov_mask;
    union_after  = ref_mask | reg_mask;
    mask_before = ref_mask & mov_mask;
    mask_after  = ref_mask & reg_mask;

    if sum(union_before(:))==0
        IoU_before = 0;
    else
        IoU_before = sum(mask_before(:)) / sum(union_before(:));
    end
    if sum(union_after(:))==0
        IoU_after = 0;
    else
        IoU_after  = sum(mask_after(:))  / sum(union_after(:));
    end
end


function g = gray_for_overlay(img, mask)
% 把任意输入（预处理灰度 / 原始 RGB）变成 overlay 用灰度：
%   组织区域亮度为 [0,1]，背景强制为 0。
% 极性由「掩膜内均值」决定，不用整图均值（整图含大片白边会误判）。

    if isempty(img)
        g = [];
        return;
    end
    if size(img,3)==3
        g = im2gray(img);
    else
        g = img;
    end
    if ~isfloat(g)
        g = im2double(g);
    elseif max(g(:)) > 1
        g = double(g) / 255;
    else
        g = double(g);
    end

    mask = mask ~= 0;
    if any(mask(:))
        if mean(g(mask)) < 0.5
            g = 1 - g;
        end
        vals = g(mask);
        lo = min(vals);
        hi = max(vals);
        if hi > lo
            g = (g - lo) / (hi - lo);
        end
    end
    g(~mask) = 0;
    g = min(max(g,0),1);
end


function overlay = compose_overlay(ref_gray, mov_gray, alpha_ref, alpha_mov)
% 红 = 参考，绿 = 移动，重叠 = 黄；背景保持黑，不会出现整幅绿底。
    overlay = zeros([size(ref_gray), 3]);
    overlay(:,:,1) = min(alpha_ref * ref_gray, 1);
    overlay(:,:,2) = min(alpha_mov * mov_gray, 1);
    overlay(:,:,3) = 0;
end
