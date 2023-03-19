function F = MovieXY(X,Y,dt,C,linecolor)
%% F = movieXY(X,Y,dt,C) �������ߣ����������ݻط�
% X ���ߺ����꣬nά������������Yͨά�ľ���
% Y ���������꣬��Ϊnά��������n*k����(k������)
% dt ���������ݵ�֮���ʱ��������λΪ��
% dtӰ�����ݻطŵĿ�����Ĭ��ֵΪ0.05s
% CΪ��ǽṹ��,�ַ����ͣ���ѡ'*' 'o'  's'��
% ����ֵFΪ�����Ŀ��frame�ļ�������������gifͼƬ
% By ZFS@wust 2021

 
hf = gcf;

% hold on

% ����ʱ������Ĭ��ֵ
if nargin == 2  || isempty(dt)
   dt = 0.05;
end

n = length(Y(:,1));
m = length(Y(1,:));
if isvector(X)        % ���X��������������չΪ��Yͬά�ľ���
    X = repmat( X,1,m );
end

if nargin < 4
    C = repmat({'*'},m,1);  % ���Ƴ�ʼ��
end

for ii = 1:m
    h(ii) = plot(X(1,ii),Y(1,ii),C{ii});     % ���Ƴ�ʼ��
end


X_1 = X(1,:);
Y_1 = Y(1,:);
dX = 0.01*( max(X) - min(X) );

Ymax = max( Y(:) );
Ymin = min( Y(:) );
Xmax = max( X(:) );
Xmin = min( X(:) );
dY = 0.05*( Ymax-Ymin );
dX = 0.05*( Xmax-Xmin );

axis([Xmin-dX Xmax+dX Ymin-dY Ymax+dY]);  % ����������

% F(1) = getframe(hf);

for ii=2:n

   for jj = 1:length(h)
      set( h(jj),'xdata',X(ii,jj),'ydata',Y(ii,jj) );  % ���µ�
   end

%    line([X_1; X(ii,:)],[Y_1;Y(ii,:)],'linewid',2,'color',[1 0 0]);   % ��������
   % �޸ģ��Ա��ò�ͬ���ߵ���ɫ�����͡��߿�
   for i =1:m
       hl = line([X_1(i); X(ii,i)],[Y_1(i);Y(ii,i)],'linewid',2,'color',linecolor(i,:));   % ��������
       set(hl,'linestyle','--');
   end
   drawnow
   X_1 = X(ii,:);
   Y_1 = Y(ii,:);
   pause(dt);
   
   F(ii) = getframe(hf);
end
