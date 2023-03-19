function  Fun_F2gif(F,filename,S)
%% Fun_F2gif(F,filename,S)  ��frame��ʽͼƬ����ļ�ת��ΪgifͼƬ
% FΪframe��ʽ�ṹ��
% SΪ���ò�����Ĭ��Ϊgif������֮֡��ļ��
% filenameΪ�ļ�����Ĭ��Ϊ'Test'
% SΪ���ò��� SΪ��֡ʱ����

% Ĭ�ϲ���
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
        imwrite(I,map,filename, 'Loopcount',inf,'DelayTime',S);      % ��֡
    else
        imwrite(I,map,filename, 'WriteMode','append','DelayTime',S); % ����֡ 
    end
    
end