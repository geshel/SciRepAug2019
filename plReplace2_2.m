load replaceAll2 E* e* nat;

% nat is in:
% million acres
% Tg or million metric tons N
% million metric tons CO2e
% million m3/year
nat = nat' .* [ 0.40469 1 1 1e-3 ]';
  
% change land units, m2 -> 1e6 ha
Emt(1,:) = Emt(1,:)*1e-10; emt(1) = emt(1)*1e-10;  vr{1} = 'land';
Ebf(1,:) = Ebf(1,:)*1e-10; ebf(1) = ebf(1)*1e-10; unt{1} = '10^6 ha'; 
eal(1)   =   eal(1)*1e-10; 

% change Nr units, g -> 1e9 kg
Emt(2,:) = Emt(2,:)*1e-12; emt(2) = emt(2)*1e-12;  vr{2} = 'nitrogen'; 
Ebf(2,:) = Ebf(2,:)*1e-12; ebf(2) = ebf(2)*1e-12; unt{2} = '10^9 kg'; 
eal(2)   =   eal(2)*1e-12; 

% change CO2eq units, g -> 1e9 kg
Emt(3,:) = Emt(3,:)*1e-12; emt(3) = emt(3)*1e-12;  vr{3} = 'emissions';
Ebf(3,:) = Ebf(3,:)*1e-12; ebf(3) = ebf(3)*1e-12; unt{3} = '10^9 kg';
eal(3)   =   eal(3)*1e-12; 

% change water units, L -> 1e9 m3
Emt(4,:) = Emt(4,:)*1e-12; emt(4) = emt(4)*1e-12;  vr{4} = 'water'; 
Ebf(4,:) = Ebf(4,:)*1e-12; ebf(4) = ebf(4)*1e-12; unt{4} = '10^9 m^3';
eal(4)   =   eal(4)*1e-12; 

% aboveUnits/person-day -> aboveUnits/nation-yr
cl  = 365*327e6; eal = eal*cl;
Emt = Emt*cl;    emt = emt*cl;
Ebf = Ebf*cl;    ebf = ebf*cl;

% resource use differences,
% useByMeat minus useByReplacementPlants
Dmt = emt*ones(1,length(Emt(1,:))) - Emt;
Dbf = ebf*ones(1,length(Ebf(1,:))) - Ebf;

% stats of resource use by the MC realizations: 
Emt = [ mean(Emt,2) std(Emt,0,2) ];
Ebf = [ mean(Ebf,2) std(Ebf,0,2) ];

% stats of resource use difference by the MC realizations: 
Dmt = [ mean(Dmt,2) std(Dmt,0,2) ];
Dbf = [ mean(Dbf,2) std(Dbf,0,2) ];

% changes as percent of national dietary resource usage 
pmtdt = round(100*(emt-Emt(:,1))./eal);
pbfdt = round(100*(ebf-Ebf(:,1))./eal);

% changes as  percent of national total resource usage 
pmttt = round(100*(emt-Emt(:,1))./nat);
pbftt = round(100*(ebf-Ebf(:,1))./nat);

% bar colors (only plant replacemernts bcz diffs)
cl      = zeros(2,3);
cl(1,:) = [  0 .7  0 ]; % beef plant replacement
cl(2,:) = [  0  1  0 ]; % all meat plant replacement

pn = 'abcd';
cr = 'Color';
un = 'Units';
fs = 'FontSize';
rt = 'Rotation';
ec = 'EdgeColor';
lw = 'LineWidth';
nr = 'normalized';
ha = 'HorizontalAlignment';

dp = [ 5 4 1 1 ];

clf

for i = 1:4
  
  axes('Position',[.09+(i-1)*.086+(i>2)*.013+(i>3)*.011 .5 .06 .32]);
  % bars
  fill(1+.47*[-1 -1 1 1],[0 Dmt(i,1)*[1 1] 0],cl(1,:),ec,cl(1,:)); hold on
  fill(2+.47*[-1 -1 1 1],[0 Dbf(i,1)*[1 1] 0],cl(2,:),ec,cl(2,:))
  % whiskers
  plot(1*[1 1],Dmt(i,1)+Dmt(i,2)*[-1 1],'k-')
  plot(2*[1 1],Dbf(i,1)+Dbf(i,2)*[-1 1],'k-')
  plot(1+.1*[-1 1],(Dmt(i,1)+Dmt(i,2))*[1 1],'k-')
  plot(1+.1*[-1 1],(Dmt(i,1)-Dmt(i,2))*[1 1],'k-')
  plot(2+.1*[-1 1],(Dbf(i,1)+Dbf(i,2))*[1 1],'k-')
  plot(2+.1*[-1 1],(Dbf(i,1)-Dbf(i,2))*[1 1],'k-')
  % axis
  mn = min([  0 (Dmt(i,1)-Dmt(i,2)) (Dbf(i,1)-Dbf(i,2)) ])*1.02;
  mx = max([ .2 (Dmt(i,1)+Dmt(i,2)) (Dbf(i,1)+Dbf(i,2)) ]);
  axis([.53 2.47 mn mx])
  % text
  if i==2 | i==3
    text(.25,.03,'all meat',rt,90,un,nr)
    text(.75,.03,  'beef'  ,rt,90,un,nr)
  end
  text(.5,1.055,unt{i},un,nr,ha,'c',fs,8.5)
  text(.5,-.05,vr{i},un,nr,ha,'c')
  text(-.1,1.06,[ '{\bf' pn(i) '}'],un,nr,ha,'c',fs,12)
  set(gca,'XTick',[],'Box','off','TickLength',[.03 .1],'YMinorTick','on')
  if i<2
    text(.4,1.15, ...
	'resource savings by replacements',un,nr)
    ylabel('\rightarrow resource savings \rightarrow')
  elseif i>3
    text(1.12,.5,'\leftarrow added resource use \leftarrow',un,nr,rt,90,ha,'c')
  end
end

cl = colormap;
cl = cl(1:12:51,:);

% percentage panel
axes('Position',[.48 .48 .4 .34])
P  = [ pbftt pbfdt pmttt pmtdt ];
ax = [ .9 16.7 min(P(:)) max(P(:)) ];
plot([1 5],[0 0],'w-'); hold on
for i = 1:4 % loop on resource
  d = P(i,:);
  for k = 1:4 % loop on above vars
    x0 = (i-1)*4.02+k*.9;
    fill(x0+.85*[0 0 1 1],[0 d(k)*[1 1] 0],cl(k,:),ec,'w')
    if d(k)>0
      for m = 5:5:d(k)
	plot(x0    +.2*[0 1],[m m],'w-')
	plot(x0+.85-.2*[0 1],[m m],'w-')
      end
    else
      for m = -5:-5:d(k)
	plot(x0    +.2*[0 1],[m m],'w-')
	plot(x0+.85-.2*[0 1],[m m],'w-')
      end
    end
  end
end
ylabel('% resource use change') 
axis(ax)
set(gca,'Box','off','YAxisLocation','r')
set(get(gca,'XAxis'),'Visible','off')
text( 2.5,-3,vr{1},ha,'c')  
text( 6.7,-3,vr{2},ha,'c')  
text(10.8,-3,vr{3},ha,'c')  
text(14.8, 3,vr{4},ha,'c')  
text(.01,.9,'{\bfe}',un,nr,ha,'c',fs,12)

% color legend for percentage panel
axes('Position',[.475 .835 .5 .04])
fill(     [0 0 3 3],[0 9 9 0],cl(1,:),ec,'w'); hold on
fill(   7+[0 0 3 3],[0 9 9 0],cl(2,:),ec,'w')
fill(14.5+[0 0 3 3],[0 9 9 0],cl(3,:),ec,'w')
fill(23.1+[0 0 3 3],[0 9 9 0],cl(4,:),ec,'w')
text(3.2,7.5,'beef,')
text(3.2,2.5,'total')
text(10.2,7.5,'beef,')
text(10.2,2.5,'dietary')
text(17.7,7.5,'all meat,')
text(17.7,2.5,'total')
text(26.3,7.5,'all meat,')
text(26.3,2.5,'dietary')
set(gca,'XLim',[0 36],'YLim',[0 9])
set(get(gca,'XAxis'),'Visible','off')
set(get(gca,'YAxis'),'Visible','off')

print -dpdf -r4350 replace2_2.pdf
