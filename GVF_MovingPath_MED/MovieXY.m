function F = MovieXY(X,Y,dt,C,linecolor)
%% F = movieXY(X,Y,dt,C) 动画曲线，可用于数据回放
% X 曲线横坐标，n维列向量，或与Y通维的矩阵
% Y 曲线纵坐标，可为n维列向量或n*k矩阵(k条曲线)
% dt 两相邻数据点之间的时间间隔，单位为秒
% dt影响数据回放的快慢，默认值为0.05s
% C为标记结构体,字符类型，可选'*' 'o'  's'等
% 返回值F为动画的框架frame文件，可用于生成gif图片
% By ZFS@wust 2021

 
hf = gcf;

% hold on

% 给出时间间隔的默认值
if nargin == 2  || isempty(dt)
   dt = 0.05;
end

n = length(Y(:,1));
m = length(Y(1,:));
if isvector(X)        % 如果X是向量，则将其扩展为与Y同维的矩阵
    X = repmat( X,1,m );
end

if nargin < 4
    C = repmat({'*'},m,1);  % 绘制初始点
end

for ii = 1:m
    h(ii) = plot(X(1,ii),Y(1,ii),C{ii});     % 绘制初始点
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

axis([Xmin-dX Xmax+dX Ymin-dY Ymax+dY]);  % 限制坐标轴

% F(1) = getframe(hf);

for ii=2:n

   for jj = 1:length(h)
      set( h(jj),'xdata',X(ii,jj),'ydata',Y(ii,jj) );  % 更新点
   end

%    line([X_1; X(ii,:)],[Y_1;Y(ii,:)],'linewid',2,'color',[1 0 0]);   % 曲线连线
   % 修改，以便获得不同曲线的颜色、线型、线宽
   for i =1:m
       hl = line([X_1(i); X(ii,i)],[Y_1(i);Y(ii,i)],'linewid',2,'color',linecolor(i,:));   % 曲线连线
       set(hl,'linestyle','--');
   end
   drawnow
   X_1 = X(ii,:);
   Y_1 = Y(ii,:);
   pause(dt);
   
   F(ii) = getframe(hf);
end
