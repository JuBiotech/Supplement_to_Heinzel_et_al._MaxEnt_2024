% plot3(solution(2,1:707),solution(3,1:707),solution(4,1:707),'s');grid on;hold on;
% plot3(solution(2,1:707),solution(3,1:707),solution(5,1:707),'o','Color','red');

% plot(solution(3,1:101),solution(4,1:101),'.-','Color','blue');
% grid on; hold on;
% plot(solution(3,1:101),solution(5,1:101),'.-','Color','red');
% 
% plot(solution(3,102:202),solution(4,102:202),'.-','Color','blue');
% plot(solution(3,102:202),solution(5,102:202),'.-','Color','green');
% 
% plot(solution(3,203:303),solution(4,203:303),'.-','Color','blue');
% plot(solution(3,203:303),solution(5,203:303),'.-','Color','red');
% 
% plot(solution(3,304:404),solution(4,304:404),'.-','Color','blue');
% plot(solution(3,304:404),solution(5,304:404),'.-','Color','green');
% 
% plot(solution(3,405:505),solution(4,405:505),'.-','Color','blue');
% plot(solution(3,405:505),solution(5,405:505),'.-','Color','red');
% 
% plot(solution(3,506:606),solution(4,506:606),'.-','Color','blue');
% plot(solution(3,506:606),solution(5,506:606),'.-','Color','green');
% 
% plot(solution(3,607:707),solution(4,607:707),'.-','Color','blue');
% plot(solution(3,607:707),solution(5,607:707),'.-','Color','red');

% colOrder = get(gca,'ColorOrder');
% 
% plot3(solution(2,1:101),solution(3,1:101),solution(4,1:101),'-', 'Color', colOrder(mod(1, size(colOrder, 1))+1,:));
% grid on; hold on;
% plot3(solution(2,1:101),solution(3,1:101),solution(5,1:101),'.', 'Color', colOrder(mod(1, size(colOrder, 1))+1,:));
% 
% plot3(solution(2,102:202),solution(3,102:202),solution(4,102:202),'-', 'Color', colOrder(mod(2, size(colOrder, 1))+1,:));
% plot3(solution(2,102:202),solution(3,102:202),solution(5,102:202),'.', 'Color', colOrder(mod(2, size(colOrder, 1))+1,:));
% 
% plot3(solution(2,203:303),solution(3,203:303),solution(4,203:303),'-', 'Color', colOrder(mod(3, size(colOrder, 1))+1,:));
% plot3(solution(2,203:303),solution(3,203:303),solution(5,203:303),'.', 'Color', colOrder(mod(3, size(colOrder, 1))+1,:));
% 
% plot3(solution(2,304:404),solution(3,304:404),solution(4,304:404),'-', 'Color', colOrder(mod(4, size(colOrder, 1))+1,:));
% plot3(solution(2,304:404),solution(3,304:404),solution(5,304:404),'.', 'Color', colOrder(mod(4, size(colOrder, 1))+1,:));
% 
% plot3(solution(2,405:505),solution(3,405:505),solution(4,405:505),'-', 'Color', colOrder(mod(5, size(colOrder, 1))+1,:));
% plot3(solution(2,405:505),solution(3,405:505),solution(5,405:505),'.', 'Color', colOrder(mod(5, size(colOrder, 1))+1,:));
% 
% plot3(solution(2,506:606),solution(3,506:606),solution(4,506:606),'-', 'Color', colOrder(mod(6, size(colOrder, 1))+1,:));
% plot3(solution(2,506:606),solution(3,506:606),solution(5,506:606),'.', 'Color', colOrder(mod(6, size(colOrder, 1))+1,:));
% 
% plot3(solution(2,607:707),solution(3,607:707),solution(4,607:707),'-', 'Color', colOrder(mod(7, size(colOrder, 1))+1,:));
% plot3(solution(2,607:707),solution(3,607:707),solution(5,607:707),'.', 'Color', colOrder(mod(7, size(colOrder, 1))+1,:));

% colOrder = get(gca,'ColorOrder');
% 
% plot(solution(3,1:101),solution(4,1:101),'-', 'Color', colOrder(mod(1, size(colOrder, 1))+1,:));
% grid on; hold on;
% plot(solution(3,1:101),solution(5,1:101),'o', 'Color', colOrder(mod(1, size(colOrder, 1))+1,:));
% 
% plot(solution(3,102:202),solution(4,102:202),'-', 'Color', colOrder(mod(2, size(colOrder, 1))+1,:));
% plot(solution(3,102:202),solution(5,102:202),'o', 'Color', colOrder(mod(2, size(colOrder, 1))+1,:));
% 
% plot(solution(3,203:303),solution(4,203:303),'-', 'Color', colOrder(mod(3, size(colOrder, 1))+1,:));
% plot(solution(3,203:303),solution(5,203:303),'o', 'Color', colOrder(mod(3, size(colOrder, 1))+1,:));
% 
% plot(solution(3,304:404),solution(4,304:404),'-', 'Color', colOrder(mod(4, size(colOrder, 1))+1,:));
% plot(solution(3,304:404),solution(5,304:404),'o', 'Color', colOrder(mod(4, size(colOrder, 1))+1,:));
% 
% plot(solution(3,405:505),solution(4,405:505),'-', 'Color', colOrder(mod(5, size(colOrder, 1))+1,:));
% plot(solution(3,405:505),solution(5,405:505),'o', 'Color', colOrder(mod(5, size(colOrder, 1))+1,:));
% 
% plot(solution(3,506:606),solution(4,506:606),'-', 'Color', colOrder(mod(6, size(colOrder, 1))+1,:));
% plot(solution(3,506:606),solution(5,506:606),'o', 'Color', colOrder(mod(6, size(colOrder, 1))+1,:));
% 
% plot(solution(3,607:707),solution(4,607:707),'-', 'Color', colOrder(mod(7, size(colOrder, 1))+1,:));
% plot(solution(3,607:707),solution(5,607:707),'o', 'Color', colOrder(mod(7, size(colOrder, 1))+1,:));

% colOrder = get(gca,'ColorOrder');
% 
% plot(solution1(3,1:101),solution1(4,1:101),'-', 'Color', 'red');
% grid on; hold on;
% plot(solution1(3,1:101),solution1(5,1:101),'.', 'Color','red');
% 
% plot(solution1(3,102:202),solution1(4,102:202),'-', 'Color','red');
% plot(solution1(3,102:202),solution1(5,102:202),'.', 'Color','red');
% 
% plot(solution1(3,203:303),solution1(4,203:303),'-', 'Color','red');
% plot(solution1(3,203:303),solution1(5,203:303),'.', 'Color','red');
% 
% plot(solution1(3,304:404),solution1(4,304:404),'-', 'Color','red');
% plot(solution1(3,304:404),solution1(5,304:404),'.', 'Color','red');
% 
% plot(solution1(3,405:505),solution1(4,405:505),'-', 'Color','red');
% plot(solution1(3,405:505),solution1(5,405:505),'.', 'Color','red');
% 
% plot(solution1(3,506:606),solution1(4,506:606),'-', 'Color','red');
% plot(solution1(3,506:606),solution1(5,506:606),'.', 'Color','red');
% 
% plot(solution1(3,607:707),solution1(4,607:707),'-', 'Color','red');
% plot(solution1(3,607:707),solution1(5,607:707),'.', 'Color','red');
% 
% plot(solution10(3,1:101),solution10(4,1:101),'-', 'Color', 'blue');
% plot(solution10(3,1:101),solution10(5,1:101),'.', 'Color','blue');
% 
% plot(solution10(3,102:202),solution10(4,102:202),'-', 'Color','blue');
% plot(solution10(3,102:202),solution10(5,102:202),'.', 'Color','blue');
% 
% plot(solution10(3,203:303),solution10(4,203:303),'-', 'Color','blue');
% plot(solution10(3,203:303),solution10(5,203:303),'.', 'Color','blue');
% 
% plot(solution10(3,304:404),solution10(4,304:404),'-', 'Color','blue');
% plot(solution10(3,304:404),solution10(5,304:404),'.', 'Color','blue');
% 
% plot(solution10(3,405:505),solution10(4,405:505),'-', 'Color','blue');
% plot(solution10(3,405:505),solution10(5,405:505),'.', 'Color','blue');
% 
% plot(solution10(3,506:606),solution10(4,506:606),'-', 'Color','blue');
% plot(solution10(3,506:606),solution10(5,506:606),'.', 'Color','blue');
% 
% plot(solution10(3,607:707),solution10(4,607:707),'-', 'Color','blue');
% plot(solution10(3,607:707),solution10(5,607:707),'.', 'Color','blue');

colOrder = get(gca,'ColorOrder');

% plot(max26Solution(3,203:303),min18Solution(4,203:303),'-', 'Color',[0.8 0.8 0.8]);
grid on; hold on;
% plot(max26Solution(3,203:303),min18Solution(5,203:303),'.', 'Color',[0.8 0.8 0.8]);
% plot(max26Solution(3,506:606),min18Solution(4,506:606),'-', 'Color',[0.8 0.8 0.8]);
% plot(max26Solution(3,506:606),min18Solution(5,506:606),'s', 'Color',[0.8 0.8 0.8]);

% plot(maxBmSolution(3,1:101),maxBmSolution(4,1:101),'-', 'Color', 'red');
% grid on; hold on;
% plot(maxBmSolution(3,1:101),maxBmSolution(5,1:101),'.', 'Color','red');

% plot(maxBmSolution(3,102:202),maxBmSolution(4,102:202),'-', 'Color','red');
% plot(maxBmSolution(3,102:202),maxBmSolution(5,102:202),'.', 'Color','red');
% 
plot(maxBmSolution(3,203:303),maxBmSolution(4,203:303),'-', 'Color','green');
plot(maxBmSolution(3,203:303),maxBmSolution(5,203:303),'.', 'Color','green');
% 
% plot(maxBmSolution(3,304:404),maxBmSolution(4,304:404),'-', 'Color','red');
% plot(maxBmSolution(3,304:404),maxBmSolution(5,304:404),'.', 'Color','red');
% 
% plot(maxBmSolution(3,405:505),maxBmSolution(4,405:505),'-', 'Color','red');
% plot(maxBmSolution(3,405:505),maxBmSolution(5,405:505),'.', 'Color','red');
% 
plot(maxBmSolution(3,506:606),maxBmSolution(4,506:606),'-', 'Color','green');
plot(maxBmSolution(3,506:606),maxBmSolution(5,506:606),'s', 'Color','green');
% 
% plot(maxBmSolution(3,607:707),maxBmSolution(4,607:707),'-', 'Color','red');
% plot(maxBmSolution(3,607:707),maxBmSolution(5,607:707),'.', 'Color','red');

% plot(max26Solution(3,1:101),max26Solution(4,1:101),'-', 'Color', 'blue');
% plot(max26Solution(3,1:101),max26Solution(5,1:101),'.', 'Color','blue');

% plot(max26Solution(3,102:202),max26Solution(4,102:202),'-', 'Color','blue');
% plot(max26Solution(3,102:202),max26Solution(5,102:202),'.', 'Color','blue');
% 
plot(max26Solution(3,203:303),max26Solution(4,203:303),'-', 'Color','black');
plot(max26Solution(3,203:303),max26Solution(5,203:303),'.', 'Color','black');

% plot(max26Solution(3,304:404),max26Solution(4,304:404),'-', 'Color','blue');
% plot(max26Solution(3,304:404),max26Solution(5,304:404),'.', 'Color','blue');
% 
% plot(max26Solution(3,405:505),max26Solution(4,405:505),'-', 'Color','blue');
% plot(max26Solution(3,405:505),max26Solution(5,405:505),'.', 'Color','blue');
% 
plot(max26Solution(3,506:606),max26Solution(4,506:606),'-', 'Color','black');
plot(max26Solution(3,506:606),max26Solution(5,506:606),'s', 'Color','black');
% 
% plot(max26Solution(3,607:707),max26Solution(4,607:707),'-', 'Color','blue');
% plot(max26Solution(3,607:707),max26Solution(5,607:707),'.', 'Color','blue');