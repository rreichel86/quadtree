clear all
elem = load('selm.txt','-ascii');
coor = load('scor.txt','-ascii');
n = size(elem,1);  % number of elements 
m = size(coor,1);  % number of nodes

xmin = min(coor(:,2));
xmax = max(coor(:,2));
ymin = min(coor(:,3));
ymax = max(coor(:,3));

f = figure('units','normalized','Resize','on','Name','QTREE','ToolBar','figure','MenuBar','none','outerposition',[0 0 0.5 1]);



for i=1:n
    
matNro = elem(i,3);    
numNodes = elem(i,2)+3;  
Nodes = elem(i,4:numNodes);
x = [coor(Nodes(1:end),2)', coor(Nodes(1),2)];
y = [coor(Nodes(1:end),3)', coor(Nodes(1),3)];

if matNro == 1
%     plot(x,y,'-r')
   patch(x,y,'red')
%     fill(x,y,'r')
elseif matNro == 2
%     plot(x,y,'-r')
    patch(x,y,'blue')
elseif matNro == 3    
%     plot(x,y,'-g')
   patch(x,y,'cyan')
else 
%    plot(x,y,'-c') 
   patch(x,y,'green')
end    


axis equal
axis ([xmin xmax ymin ymax])
axis off

%pause(0.1)
hold on
end 


plot(coor(m-n+1:m,2),coor(m-n+1:m,3),'k*')
