function y = sig(u,m )
% SIG y = sig(u,m )
% date: 6th March 2022
% Author: quyinsong
% Reference: 固定时间预测器下的欠驱动无人艇路径跟踪控制

y = abs(u)^m*sign(u);

end

