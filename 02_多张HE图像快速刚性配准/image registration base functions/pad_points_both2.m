function pts_out = pad_points_both2(pts, im_size, target_size, ext)
    % 与 pad_im_both2 等价的“点坐标版”padding
    % pts: N×2，[x,y] 像素坐标（列、行）
    % im_size:   [H W] 原图像的尺寸（pad/crop 之前）
    % target_size: [H W] 目标尺寸（= refsize）
    % ext: 标量，等价于 pad_im_both2 的 ext（= padall）
    
        if nargin < 4, ext = 0; end
    
        % 若原图比 target 大，配准代码会先裁剪到 min(im_size, target_size)
        base_size = [min(im_size(1), target_size(1)), ...
                     min(im_size(2), target_size(2))];
    
        % 计算需要补齐到 target 的总 padding 量，并按 pad_im_both2 的规则拆到两边
        szim = [target_size(1) - base_size(1), ...
                target_size(2) - base_size(2)];
        szA = floor(szim/2);        % pre（上/左）
        szB = szim - szA + ext;     % post（下/右）— 这里与图像一致，但用不着数值
        szA = szA + ext;            % pre 再加 ext（注意：两边都加 ext）
    
        % 点坐标加“前侧”(pre) padding 量
        % 注意：x 加列方向的 padding（szA(2)），y 加行方向的 padding（szA(1)）
        pts_out = pts + [szA(2), szA(1)];
    end
    