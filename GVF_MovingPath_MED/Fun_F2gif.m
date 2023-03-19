function  Fun_F2gif(F,filename,S)
%% Fun_F2gif(F,filename,S)  将frame格式图片框架文件转换为gif图片
% F为frame格式结构体
% S为设置参数，默认为gif动画两帧之间的间隔
% filename为文件名，默认为'Test'
% S为设置参数 S为单帧时间间隔

% 默认参数
if nargin < 2
    filename = 'Test.gif';
    S = 0.1;
elseif nargin < 3
    S = 0.1;
elseif nargin == 3
    if isempty(filename)
        filename = 'Test.gif';
    end
end

for ii = 1:length(F)
    if iscell(F)
       f = F{ii};
       [I,map] = rgb2ind(f,256);
    else
       f = F(ii); 
       I = frame2im(f);
       [I,map] = rgb2ind(I,256);
    end

    if ii == 1
        imwrite(I,map,filename, 'Loopcount',inf,'DelayTime',S);      % 首帧
    else
        imwrite(I,map,filename, 'WriteMode','append','DelayTime',S); % 后续帧 
    end
    
end