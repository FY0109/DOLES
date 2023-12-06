x=[1 2 3 4];
 y=[log(664) log(654) log(412) log(20) log(349) log(112) log(135);
     log(2058) log(1479) log(1310) log(67) log(1097) log(177) log(189);
     log(6420) log(8970) log(13590) log(176) log(2350) log(450) log(506);
     log(12882) log(14562) log(26788) log(316) log(3350) log(939) log(1031)]; 
 figureHandle = figure;
 bar(x,y);
 hYLabel1 = ylabel('Logarithm of Running Time ');
 hLegend = legend( 'FPMVS','MSGL' ,'SMVSC' ,'OPMC', 'OMSC','AWMVC','Ours');

set(gca, 'FontName', 'Helvetica', 'FontSize', 10)
 
set(gca,'xTicklabel',{'YTF10','YTF20','YTF50','YTF100'})
fileout = 'test';
print(figureHandle,[fileout,'.png'],'-r600','-dpng');
