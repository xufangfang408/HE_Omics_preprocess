function imG=register_global_im(im,tform,cent,f,fillval)
% applies previously calculated global image registration to an image
% Written in 2020 by Ashley Lynn Kiemen, Johns Hopkins University
% please cite Kiemen et al, Nature methods (2022)
% last updated in December 2023 by ALK

    Rin=imref2d(size(im));
    Rin.XWorldLimits = Rin.XWorldLimits-cent(1);
    Rin.YWorldLimits = Rin.YWorldLimits-cent(2);

    if f==1
        im=im(end:-1:1,:,:);
    end

    % RGB 的 FillValues 必须是长度为 3 的向量，不能用 1x1x3。
    % 1x1x3 会触发 imwarp: size must match dimensions 3 to N of A。
    fv=format_fillvalues(fillval, size(im,3));
    imG=imwarp(im,Rin,tform,'nearest','Outputview',Rin,'Fillvalues',fv);
end


function fv=format_fillvalues(fillval, nchan)
    if nargin<1 || isempty(fillval)
        if nchan==3
            fv=[0 0 0];
        else
            fv=0;
        end
        return;
    end
    v=double(fillval(:));
    if nchan==3
        if numel(v)==1
            v=repmat(v,3,1);
        elseif numel(v)<3
            v(end+1:3,1)=v(end);
        end
        fv=v(1:3)';  % 1x3，与 size(A,3)=3 匹配
    else
        fv=v(1);
    end
end
