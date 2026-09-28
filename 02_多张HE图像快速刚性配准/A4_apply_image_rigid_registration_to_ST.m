function A4_apply_image_rigid_registration_to_ST(ST_position_pth,pthdata,scale)

% 2025/12/10 更新：删掉所有弹性配准相关的代码，并优化代码结构，使其更加简洁

    % 设置路径
    pthim = "/disk/sdc/xufangfang/workdir/CODA_tif/dataset/Schmidtea_mediterranea/registered_data/24h/HE/registeredE"; % 配准后图像路径
    pthdata = "/disk/sdc/xufangfang/workdir/CODA_tif/dataset/Schmidtea_mediterranea/registered_data/24h/HE/downsample_image/registered/elastic registration/save_warps";
    pthposition = "/disk/sdc/xufangfang/workdir/CODA_tif/dataset/Schmidtea_mediterranea/registered_data/24h/locs";
    scale = 2;  % 假设缩放比例为 3，可根据实际调整

    pthim = char(pthim);
    pthdata = char(pthdata);
    pthposition = char(pthposition);


    % 获取 TSV 文件列表
    tsvlist = dir([pthposition, '/*.tsv']);

    % 加载参考尺寸和填充参数
    matlist = dir([pthdata, '/*.mat']);
    try
        datafileE = [pthdata, '/', matlist(1).name];
        load(datafileE, 'szz', 'padall');
    catch
        datafileE = [pthdata, '/', matlist(end).name];
        load(datafileE, 'szz', 'padall');
    end

    padall = ceil(padall * scale);  % 填充量缩放
    refsize = ceil(szz * scale);    % 参考图像尺寸缩放

    % 加载裁剪和旋转参数（如果存在）
    outpth = pthim;
    if exist([outpth, 'crop_data.mat'], 'file')
        load([outpth, 'crop_data.mat'], 'rot', 'rr');
    else
        rot = 0;  % 无旋转
        rr = [1, 1, refsize(2), refsize(1)];  % 无裁剪
    end

    % 遍历每个 TSV 文件
    for kz = 1:length(tsvlist)
        tsvname = tsvlist(kz).name;
        matname = [tsvname(1:end-4), '.mat'];  % TSV 和 MAT 文件名一致
        imnm = [tsvname(1:end-4), '.jpg'];
        datafileE = [pthdata, '/', matname];
        
        dataE = load(datafileE, 'tform', 'cent', 'f');

        % 如果不是参考图像，直接配准
        if isfield(dataE, 'tform') && isfield(dataE, 'cent') && isfield(dataE, 'f')

            % 加载配准参数
            load(datafileE, 'tform', 'cent', 'f');

            % 加载 TSV 文件中的点坐标
            data = readtable([pthposition, '/', tsvname], 'FileType', 'text', 'Delimiter', '\t');
            points = [data.HE_X, data.HE_Y];
            
            % —— 使用与 pad_im_both2 等价的点 padding —— 
            % 注意：需要原图（未 pad 未注册版本）的尺寸来计算“前侧/后侧”各自 padding 量
            % 原始（未注册、未 pad）的 H&E 路径 = registeredE 的上一级目录
            pthim_raw = fileparts(pthim);
            imnm_noext = tsvname(1:end-4);
            im_raw = [];
            try
                im_raw = imread(fullfile(pthim_raw, [imnm_noext '.jpg']));
            catch
            end
            if isempty(im_raw)
                try
                    im_raw = imread(fullfile(pthim_raw, [imnm_noext '.tif']));
                catch
                end
            end
            if isempty(im_raw)
                try
                    im_raw = imread(fullfile(pthim_raw, [imnm_noext '.jp2']));
                catch
                end
            end
            assert(~isempty(im_raw), '找不到原始 H&E 图像（未注册版本），无法为点计算与图像一致的 padding。');

            orig_size = size(im_raw(:,:,1));           % [H W]
            points = pad_points_both2(points, orig_size, refsize, padall);

            % 步骤 1：预处理 - 填充和翻转
            if f == 1
                points(:,2) = refsize(1) - points(:,2) + 1;  % 上下翻转
            end
            
            % 步骤 2：应用仿射变换 tform
            cent = cent * scale;  % 缩放中心点
            tform.T(3,1:2) = tform.T(3,1:2) * scale;  % 缩放平移分量
            % 调整坐标系到变换中心
            points(:,1) = points(:,1) - cent(1);
            points(:,2) = points(:,2) - cent(2);
            % 应用变换
            points_transformed = transformPointsForward(tform, points);
            % 恢复坐标系
            points_transformed(:,1) = points_transformed(:,1) + cent(1);
            points_transformed(:,2) = points_transformed(:,2) + cent(2);

            % 步骤 4：应用旋转和裁剪
            if rot ~= 0
                % 旋转
                theta = rot;
                R = [cosd(theta), -sind(theta); sind(theta), cosd(theta)];
                center = [refsize(2)/2, refsize(1)/2];
                points_transformed = (points_transformed - center) * R' + center;
            end
            if ~isempty(rr)
                % 裁剪
                points_transformed(:,1) = points_transformed(:,1) - rr(1) + 1;
                points_transformed(:,2) = points_transformed(:,2) - rr(2) + 1;
            end

            % 步骤 5：保存变换后的坐标
            data.HE_X_transformed = points_transformed(:,1);
            data.HE_Y_transformed = points_transformed(:,2);

        else
            % 参考图像没有刚性配准的参数
            data = readtable([pthposition, '/', tsvname], 'FileType', 'text', 'Delimiter', '\t');
            points = [data.HE_X, data.HE_Y];

            % —— 使用与 pad_im_both2 等价的点 padding —— 
            % 注意：需要原图（未 pad 未注册版本）的尺寸来计算“前侧/后侧”各自 padding 量
            % 原始（未注册、未 pad）的 H&E 路径 = registeredE 的上一级目录
            pthim_raw = fileparts(pthim);
            imnm_noext = tsvname(1:end-4);
            im_raw = [];
            try
                im_raw = imread(fullfile(pthim_raw, [imnm_noext '.jpg']));
            catch
            end
            if isempty(im_raw)
                try
                    im_raw = imread(fullfile(pthim_raw, [imnm_noext '.tif']));
                catch
                end
            end
            if isempty(im_raw)
                try
                    im_raw = imread(fullfile(pthim_raw, [imnm_noext '.jp2']));
                catch
                end
            end
            assert(~isempty(im_raw), '找不到原始 H&E 图像（未注册版本），无法为点计算与图像一致的 padding。');

            orig_size = size(im_raw(:,:,1));           % [H W]
            points = pad_points_both2(points, orig_size, refsize, padall);

            data.HE_X_transformed = points(:, 1);
            data.HE_Y_transformed = points(:, 2);
        end
        output_file = [pthposition, '/', tsvname];
        writetable(data, output_file, 'FileType', 'text', 'Delimiter', '\t');
        fprintf('已处理并保存：%s\n', output_file);
    end

end
