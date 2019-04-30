% Model G - Crank Nicolson Solver - plotting output from c++ solver file : modelg1d_bjd4d.cpp

% Load .txt file arrays of time and etheron potentials as 2d arrays.
%{
A = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_G.txt')
B = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_X.txt')
C = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_Y.txt')
D = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_lapG.txt')
E = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_lapX.txt')
F = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_lapY.txt')
%}
G = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_G_r.txt')
H = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_X_r.txt')
I = load('/home/bjd/software/SQK_2/crank_nicolson/model_g_1d-master/modelg_1d_outfile_Y_r.txt')



%==========================================================

figure('Position', [100, 100, 1049, 910]);

% Assign array to radius r and potentials: G, X, Y
r=G(:,1); G=G(:,2); X=H(:,2); Y=I(:,2);

%producing spline fit to data
xq2 = 0:1:800;

s7 = spline(r,G,xq2);
s8 = spline(r,X,xq2);
s9 = spline(r,Y,xq2);

%plot(t,G,'o',xq1,s,'-.',t,X,'x',xq1,s2,'-',t,Y,'s',xq1,s3,'.')

%plot(t,G,'o',xq1,s,'-.',t,X,'x',xq1,s2,'-',t,Y,'s',xq1,s3,'.')
%plot(t,G,xq1,s,'-.',t,X,xq1,s2,'-',t,Y,xq1,s3,'.')

%plot(t,G,xq1,s,'-.',t,X,xq1,s2,'-',t,Y,xq1,s3,'.',t,lapG,xq1,s4,'o',t,lapX,xq1,s5,'s',t,lapY,xq1,s6,'x')
plot(r,G,xq2,s7,'.','color','green',r,X,xq2,s8,'.','color','red',r,Y,xq2,s9,'.','color','black')

%scatter (r, G);

%title ("scatter() plot with red bubbles");%legend('Sample Points','pchip','spline','Location','SouthEast')
legend('G','X','Y','Location','East')

title({'Model G - Crank Nicolson Solver';'radius vs. potentials: G, X, Y- sline plot'})

xlabel('radius / dimensionless')
ylabel('G, X, Y / concentration / dimensionless')

%saveas(gcf,'/home/bjd/octave/examples/code/bjd/Model G - Crank Nicolson Solver 1g.png')
saveas(1,'/home/bjd/octave/examples/code/bjd/Model G - Crank Nicolson Solver radius vs potentials 1a.png')
%==========================================================

%==========================================================

%==========================================================

%==========================================================

