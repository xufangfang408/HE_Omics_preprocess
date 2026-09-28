function A3_apply_image_rigid_registration(pthim_in,pthdata_in,scale_in,padnum,cropim,redo,outpth_in)

% 2025/12/10 更新：删掉所有弹性配准相关的代码
% 2026/9/11  更新：填充色改用背景中位数（原来的 nested mode 会把深染切片填成洋红）；
%                  修正 base functions 的 addpath；输出目录提到开头统一配置；
%                  路径参数改为「开头默认值 + 入参可覆盖」，与 A2 一致

    pthim = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/interim/PGD_dataset/HE_unregistered";
    pthdata = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/interim/PGD_dataset/HE_unregistered/downsample_image/registered/save_warps";
    scale = 2;
    % 输出目录；置为 "" 则退回旧行为 [pthim,'registeredE/']
    % outpth = "/mnt/zzf_nas/A_xff/Spaceland-omics/data/processed/PGD_dataset/HE_registered";
    outpth = "";

    if nargin>=1 && ~isempty(pthim_in); pthim=pthim_in; end
    if nargin>=2 && ~isempty(pthdata_in); pthdata=pthdata_in; end
    if nargin>=3 && ~isempty(scale_in); scale=scale_in; end
    if nargin>=7 && ~isempty(outpth_in); outpth=outpth_in; end

    pthim = char(pthim);  % 确保是 char 类型
    pthdata = char(pthdata);
    outpth = char(outpth);

    if ~exist('redo','var');redo=0;end
    if ~exist('padnum','var');pd=1;padnum=[];else;pd=0;end
    if isempty(padnum);pd=1;end
    if ~exist('cropim','var');cropim=0;end

    % add base functions to the MATLAB search path
    path(path,'/mnt/zzf_nas/A_xff/Spaceland-omics/spaceland_omics/rigid_registration/image registration base functions');

    if pthim(end)~='/';pthim=[pthim,'/'];end
    if pthdata(end)~='/';pthdata=[pthdata,'/'];end
    imlist=dir([pthim,'*tif']);fl='tif';
    if isempty(imlist);imlist=dir([pthim,'*jp2']);fl='jp2';end
    if isempty(imlist);imlist=dir([pthim,'*jpg']);fl='jpg';end
    if isempty(outpth);outpth=[pthim,'registeredE/'];end
    if outpth(end)~='/';outpth=[outpth,'/'];end
    if ~isfolder(outpth);mkdir(outpth);end
    matlist=dir([pthdata,'*mat']);

    try 
        datafileE=[pthdata,matlist(1).name];
        load(datafileE,'szz','padall');
    catch
        datafileE=[pthdata,matlist(end).name];
        load(datafileE,'szz','padall');
    end

    padall=ceil(padall*scale);
    refsize=ceil(szz*scale);

    % determine crop region
    if cropim~=0
        if exist([outpth,'crop_data.mat'],'file')
            load([outpth,'crop_data.mat'],'rot','rr');
        else
            if length(cropim)==1
                [rot,rr]=get_cropim(pthdata,scale);
            else
                rot=cropim(1);rr=cropim(2:end);
            end
            save([outpth,'crop_data.mat'],'rot','rr');
        end
    end

    % register each image and save to outpth
    count=1;
    for kz=1:length(matlist)
        imnm=[matlist(kz).name(1:end-3),fl];outnm=imnm;
        disp(['registering image ',num2str(kz),' of ',num2str(length(matlist)),': ',imnm])
        if exist([outpth,outnm],'file') && ~redo;disp('  already registered');continue;end
        
        
        if ~exist([pthim,imnm],'file');continue;end
        datafileE=[pthdata,imnm(1:end-3),'mat'];
        
        % load image
        IM=imread([pthim,imnm]);
        szim=size(IM(:,:,1));
        if pd;padnum=background_fillval(IM);end
        if szim(1)>refsize(1) || szim(2)>refsize(2)
            a=min([szim; refsize]);
            IM=IM(1:a(1),1:a(2),:);
        end
        IM=pad_im_both2(IM,refsize,padall,padnum);
        
        % if not reference image, register
        try
            load(datafileE,'tform','cent','f');
            if f==1;IM=IM(end:-1:1,:,:);end
            
            % === tform不应该包含缩放，进行断言判断 ===
            % 提取线性部分
            A = tform.T(1:2,1:2);
            % 计算条件
            col1_norm = norm(A(:,1));
            col2_norm = norm(A(:,2));
            orthogonality = dot(A(:,1), A(:,2));
            detA = det(A);
            % 设置容差
            tol = 1e-3;
            % 断言：列长度=1，正交，行列式=1
            assert(abs(col1_norm-1) < tol && ...
                abs(col2_norm-1) < tol && ...
                abs(orthogonality) < tol && ...
                abs(detA-1) < tol, ...
                'Error: tform 包含缩放或非刚性变换！');
            % ==============================

            IMG=register_IM(IM,tform,scale,cent,padnum);
            
            % 如果不使用弹性配准
            IME = IMG;
            
        catch
            IME=IM;
        end
        
        
        if cropim
            IME=imrotate(IME,rot,'nearest');
            IME=imcrop(IME,rr);
        end
        
        imwrite(IME,[outpth,outnm]);

        count=count+1;
        disp('  done');
        clearvars tform rsc cent D f
    end

end

function fillval=background_fillval(im)
% 填充色 = 图像四周边框带的逐通道中位数，与 preprocessing.m 中同名函数的无掩膜分支一致。
%
% 不能用 squeeze(mode(mode(im,2),1))'：该写法逐通道独立取众数，而深染切片的绿通道
% 常被裁剪到 1（S050_2 有 14% 的像素 G==1，比任何一个背景绿值都多），于是逐行众数
% 是 1，整体众数也是 1，填充色变成 [156 1 154] 这类洋红。中位数不受该尖峰影响。
    if size(im,3)==1;im=repmat(im,[1 1 3]);end
    [h,w,~]=size(im);
    bw=max(8,round(min(h,w)*0.04));
    fillval=zeros(1,3);
    for c=1:3
        ch=im(:,:,c);
        band=[reshape(ch(1:bw,:),[],1);reshape(ch(end-bw+1:end,:),[],1);...
              reshape(ch(:,1:bw),[],1);reshape(ch(:,end-bw+1:end),[],1)];
        fillval(c)=median(double(band));
    end
end


function IM=register_IM(IM,tform,scale,cent,abc)    
    % rough registration
    cent=cent*scale;
    tform.T(3,1:2)=tform.T(3,1:2)*scale;
    Rin=imref2d(size(IM));
        Rin.XWorldLimits = Rin.XWorldLimits-cent(1);
        Rin.YWorldLimits = Rin.YWorldLimits-cent(2);
    IM=imwarp(IM,Rin,tform,'nearest','outputview',Rin,'fillvalues',abc);
end


function [rot,rr]=get_cropim(pthdata,scale)

    pth1=[pthdata,'../'];
    imlist=dir([pth1,'*tif']);if isempty(imlist);imlist=dir([pth1,'*jpg']);end
    im1=rgb2gray(imread([pth1,imlist(1).name]));
    im2=rgb2gray(imread([pth1,imlist(round(length(imlist)/2)).name]));
    im3=rgb2gray(imread([pth1,imlist(end).name]));
    im=cat(3,im1,im2,im3);
    h=figure;imshow(im);isgood=0;
    while isgood~=1
        rot=input('angle?\n');
        imshow(imrotate(im,rot));
        isgood=input('is good?\n');
    end
    im=imrotate(im,rot,'nearest');
    [~,rr]=imcrop(im);
    rr=round(rr)*scale;
    close(h)

end


