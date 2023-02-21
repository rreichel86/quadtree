clear all
close all
elem = load('./Output/plot/selm.txt','-ascii');
coor = load('./Output/plot/scor.txt','-ascii');
n = size(elem,1);  % number of elements 
m = size(coor,1);  % number of nodes

xmin = min(coor(:,2));
xmax = max(coor(:,2));
ymin = min(coor(:,3));
ymax = max(coor(:,3));

f = figure('units','normalized','Resize','on','Name','QTREE MESH','ToolBar','figure','MenuBar','none','outerposition',[0 0 0.5 1]);

axis equal
axis ([xmin xmax ymin ymax])
axis off

hold on

matNro = zeros(1,n);
numNodes = zeros(1,n);

matNro = elem(:,3);
maxNumMat = max(matNro);
numNodes = elem(:,2);
maxNumNodes = max(numNodes);

elemNodes = zeros(n,maxNumNodes);
elemNodes = elem(:,4:3+maxNumNodes);

for i = 1:n
    patch('Faces',elemNodes(i,1:numNodes(i)),'Vertices',coor(:,2:3),...
          'FaceVertexCData',matNro(i),'FaceColor','flat');
end
hold off



fileName = input("Please enter file name: [QtreeMesh]: ", "s");

if isempty(fileName)
    fileName = 'QtreeMesh';
    print(fileName,'-dpng');
else 
    print(fileName,'-dpng');
end



% plot(coor(m-n+1:m,2),coor(m-n+1:m,3),'k*')
