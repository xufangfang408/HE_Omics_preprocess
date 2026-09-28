function [he_ds_dir, mask_ds_dir] = downsample_he_and_masks(he_dir, mask_dir, ds_name, factor, redo)
% 将 HE 与对应 mask 用同一目标尺寸、最近邻降采样，保证组织轮廓对齐。
%
% 输入 mask 查找顺序：mask_dir/<fname>，然后 mask_dir/mask/<fname>
% 输出：
%   he_ds_dir   = HE_dir/<ds_name>/
%   mask_ds_dir = mask_dir/<ds_name>/
% 同时把二值 mask 写到 he_ds_dir/TA/<stem>.tif，供 get_ims 直接读取。

    if nargin < 5 || isempty(redo); redo = false; end
    if isempty(factor) || factor < 1
        error('downsample_factor 必须 >= 1');
    end

    he_dir = ensure_slash(char(he_dir));
    mask_dir = ensure_slash(char(mask_dir));
    ds_name = char(ds_name);
    if isempty(ds_name)
        error('downsample_name 不能为空');
    end

    he_ds_dir = ensure_slash([he_dir, ds_name]);
    mask_ds_dir = ensure_slash([mask_dir, ds_name]);
    ta_dir = ensure_slash([he_ds_dir, 'TA']);
    if ~isfolder(he_ds_dir); mkdir(he_ds_dir); end
    if ~isfolder(mask_ds_dir); mkdir(mask_ds_dir); end
    if ~isfolder(ta_dir); mkdir(ta_dir); end

    tifs = dir([he_dir, '*.tif']);
    if isempty(tifs)
        tifs = dir([he_dir, '*.tiff']);
    end
    if isempty(tifs)
        error('HE_dir 中没有 .tif: %s', he_dir);
    end

    fprintf('Downsampling HE + mask  (factor=%g)  %d files\n', factor, numel(tifs));
    fprintf('  HE   : %s -> %s\n', he_dir, he_ds_dir);
    fprintf('  mask : %s -> %s\n', mask_dir, mask_ds_dir);

    for k = 1:numel(tifs)
        fname = tifs(k).name;
        he_in = [he_dir, fname];
        he_out = [he_ds_dir, fname];
        mask_in = resolve_mask_file(mask_dir, fname);
        mask_out = [mask_ds_dir, fname];
        ta_out = [ta_dir, strip_ext(fname), '.tif'];

        already = exist(he_out, 'file') && exist(mask_out, 'file') && exist(ta_out, 'file');
        if already && ~redo
            fprintf('  [%d/%d] skip %s (exists)\n', k, numel(tifs), fname);
            continue;
        end

        tic;
        fprintf('  [%d/%d] %s\n', k, numel(tifs), fname);
        img = imread(he_in);
        if size(img, 3) > 3
            img = img(:, :, 1:3);
        elseif size(img, 3) == 1
            img = cat(3, img, img, img);
        end

        he_hw = size(img(:, :, 1));
        new_size = max(1, floor(he_hw / factor));

        TA = read_binary_mask(mask_in);
        if ~isequal(size(TA), he_hw)
            error('mask 与 HE 尺寸不一致: %s  HE=[%d %d]  mask=[%d %d]', ...
                fname, he_hw(1), he_hw(2), size(TA, 1), size(TA, 2));
        end

        % 目标尺寸由 HE 决定，mask 强制 resize 到同一 new_size，避免各自 floor 后差 1 像素
        img_ds = imresize(img, new_size, 'nearest');
        ta_ds = imresize(uint8(TA) * 255, new_size, 'nearest') > 127;

        imwrite(img_ds, he_out);
        imwrite(uint8(ta_ds) * 255, mask_out);
        imwrite(uint8(ta_ds) * 255, ta_out);
        fprintf('    %dx%d -> %dx%d  tissue=%.1f%%  %.1fs\n', ...
            he_hw(1), he_hw(2), new_size(1), new_size(2), 100 * mean(ta_ds(:)), toc);
    end
    fprintf('Downsampling complete.\n');
end


function p = ensure_slash(p)
    p = char(p);
    if isempty(p)
        error('路径为空');
    end
    if p(end) ~= '/' && p(end) ~= '\'
        p = [p, '/'];
    end
end


function stem = strip_ext(fname)
    [~, stem, ~] = fileparts(fname);
end


function p = resolve_mask_file(mask_dir, fname)
    mask_dir = ensure_slash(mask_dir);
    cands = {
        [mask_dir, fname]
        [mask_dir, 'mask/', fname]
        };
    for i = 1:numel(cands)
        if exist(cands{i}, 'file')
            p = cands{i};
            return;
        end
    end
    error('找不到与 %s 对应的 mask。已试:\n  %s\n  %s', fname, cands{1}, cands{2});
end


function TA = read_binary_mask(path)
    TA = imread(path);
    if ndims(TA) == 3
        TA = TA(:, :, 1);
    end
    TA = TA > 0;
end
