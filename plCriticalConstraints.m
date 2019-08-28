load replaceAll2 C D X* w sineq Ipp Nst

% combined normalized costs, [m2, gNr, gCO2e, L]/g
cc =  (w(1)*(C(:,1)'-mean(C(:,1)))./std(C(:,1))) ...
    + (w(2)*(C(:,2)'-mean(C(:,2)))./std(C(:,2))) ...  
    + (w(3)*(C(:,3)'-mean(C(:,3)))./std(C(:,3))) ...
    + (w(4)*(C(:,4)'-mean(C(:,4)))./std(C(:,4)));
cc  = cc-min(cc);

xbf = mean(Xbf,2);
xmt = mean(Xmt,2);

% keep only constraints critical+non-redundant to beef
s   = std(D,[],2);
i   = find( s>0 & Nst(:,2)<.1 );
D   = D(i,:); 
snq = sineq(i);    clear X* vrpp sineq Nst

% normalized diff of each plant item in each constraint 
i   = ones(size(D(1,:)));
Dn  = ( D - mean(D,2)*i ) ./ (std(D,0,2)*i);

% upper bound constraints must be negative
i       = find(snq==1);
Dn(i,:) = -Dn(i,:);

% sum over all critical constraints
dn = sum(Dn)';

[xbf,i] = sort(xbf,'descend'); 
I   = Ipp(i);
dn  = dn(i)-min(dn);
cc  = cc(i);
 
disp('------- beef top 10 items ---------')
i = 1:10;
disp([' mean daily mass = ' num2str(mean(xbf(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- beef items 10-20 ---------')
i = 11:20;
disp([' mean daily mass = ' num2str(mean(xbf(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- beef items 21-64 ---------')
i = 21:64;
disp([' mean daily mass = ' num2str(mean(xbf(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- beef items 21-54 ---------')
i = 21:54;
disp([' mean daily mass = ' num2str(mean(xbf(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- beef items 55-64 ---------')
i = 55:64;
disp([' mean daily mass = ' num2str(mean(xbf(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])
disp('---------------------------------------------------')
disp('---------------------------------------------------')

Cm  = brewermap(64,'*Spectral');  colormap(Cm)
x   = linspace(min(dn),max(dn),64)';
Cl  = interp1(x,Cm,dn);

clf

axes('Position',[.07 .3 .83 .61])
plot3(1:length(cc),xbf,cc*0,'k-','LineWidth',.1); hold on
for i = 1:10:40
  for j = 0:9
    k = i+j;
    text(i+(i-1)*.3,30,3.3-j*.2,[ int2str(k) '. ' I{k}], ...
	'Rotation',11,'FontSize',7)
  end
end
for j = 50:59
  text(length(cc)+.9,28,3.3-(j-50)*.2,[ int2str(j) '. ' I{j}], ...
      'Rotation',-17.5,'FontSize',7)
end
for j = 60:length(cc)
  text(length(cc)+.9,18,3.3-(j-60)*.2,[ int2str(j) '. ' I{j}], ...
      'Rotation',-17.5,'FontSize',7)
end
for i = 1:length(cc)
  % plot projections of the bars onto the y = 30 plane
  fill3(i+.3*[-1 -1 1 1],30*[1 1 1 1],[0 cc(i)*[1 1] 0], ...
      Cl(i,:),'EdgeColor',Cl(i,:),'FaceAlpha',.5,'EdgeAlpha',.5)
  % plot the color bars themselves
  h = plot3([i i],xbf(i)*[1 1],[0 cc(i)],'-','LineWidth',3);
  set(h,'Color',Cl(i,:))
  % plot environ cost white tickmarks on the bars
  for j = .1:.1:cc(i)
    plot3(i+.4*[-1.1 .6],xbf(i)-.04*[1 1]-.05,[j j],'w-','LineWidth',.1)
  end
end
grid
axis([0 length(cc)+1 -.5 30 0 3.5])
text(0,-8,0,'plant items arrranged by mass in diet','Rotation',11)
text(-13,21,0,'mass in replacement diet, g d^{-1}','Rotation',-20)
zlabel('combined environmental cost')
set(gca,'XTick',0:10:70,'YTick',0:5:30,'ZTick',0:.5:4)

axes('Position',[.91 .4 .02 .5])
contourf([1 2],x,[x x],x,'w')
set(gca,'YAxisLocation','r','XTick',[],'YTick',0:5:90,'Box','off')
ylabel('relative suitability for satisfyiing constraints')

print -dpdf -r350 suitability.pdf

%%%%%%%

clf

% combined normalized costs, [m2, gNr, gCO2e, L]/g
cc =  (w(1)*(C(:,1)'-mean(C(:,1)))./std(C(:,1))) ...
    + (w(2)*(C(:,2)'-mean(C(:,2)))./std(C(:,2))) ...  
    + (w(3)*(C(:,3)'-mean(C(:,3)))./std(C(:,3))) ...
    + (w(4)*(C(:,4)'-mean(C(:,4)))./std(C(:,4)));
cc  = cc-min(cc);

% sum over all critical constraints
dn = sum(Dn)';

[xmt,i] = sort(xmt,'descend'); 
I   = Ipp(i);
dn  = dn(i)-min(dn);
cc  = cc(i);
 
disp('------- all meat top 10 items ---------')
i = 1:10;
disp([' mean daily mass = ' num2str(mean(xmt(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- all meat items 10-20 ---------')
i = 11:20;
disp([' mean daily mass = ' num2str(mean(xmt(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- all meat items 21-end ---------')
i = 21:64;
disp([' mean daily mass = ' num2str(mean(xmt(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- all meat items 21-54 ---------')
i = 21:54;
disp([' mean daily mass = ' num2str(mean(xmt(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])

disp('------- all meat items 55-end ---------')
i = 55:64;
disp([' mean daily mass = ' num2str(mean(xmt(i))) ])
disp([' mean environ cost = ' num2str(mean(cc(i))) ])
disp([' mean suitability index = ' num2str(mean(dn(i))) ])
disp('---------------------------------------------------')
disp('---------------------------------------------------')

Cm  = brewermap(64,'*Spectral');  colormap(Cm)
x   = linspace(min(dn),max(dn),64)';
Cl  = interp1(x,Cm,dn);

clf

axes('Position',[.07 .3 .83 .61])
plot3(1:length(cc),xmt,cc*0,'k-','LineWidth',.1); hold on
for i = 1:10:40
  for j = 0:9
    k = i+j;
    text(i+(i-1)*.3,38,3.3-j*.2,[ int2str(k) '. ' I{k}], ...
	'Rotation',11,'FontSize',7)
  end
end
for j = 50:59
  text(length(cc)+.9,28,3.3-(j-50)*.2,[ int2str(j) '. ' I{j}], ...
      'Rotation',-17.5,'FontSize',7)
end
for j = 60:length(cc)
  text(length(cc)+.9,18,3.3-(j-60)*.2,[ int2str(j) '. ' I{j}], ...
      'Rotation',-17.5,'FontSize',7)
end
for i = 1:length(cc)
  % plot projections of the bars onto the y = 30 plane
  fill3(i+.3*[-1 -1 1 1],38*[1 1 1 1],[0 cc(i)*[1 1] 0], ...
      Cl(i,:),'EdgeColor',Cl(i,:),'FaceAlpha',.5,'EdgeAlpha',.5)
  % plot the color bars themselves
  h = plot3([i i],xmt(i)*[1 1],[0 cc(i)],'-','LineWidth',3);
  set(h,'Color',Cl(i,:))
  % plot environ cost white tickmarks on the bars
  for j = .1:.1:cc(i)
    plot3(i+.4*[-1.1 .6],xmt(i)-.04*[1 1]-.05,[j j],'w-','LineWidth',.1)
  end
end
grid
axis([0 length(cc)+1 -.5 38 0 3.5])
text(2,-10,0,'plant items arrranged by mass in diet','Rotation',11)
text(-13,21,0,'mass in replacement diet, g d^{-1}','Rotation',-20)
zlabel('combined environmental cost')
set(gca,'XTick',0:10:70,'YTick',0:5:30,'ZTick',0:.5:4)

axes('Position',[.91 .4 .02 .5])
contourf([1 2],x,[x x],x,'w')
set(gca,'YAxisLocation','r','XTick',[],'YTick',0:5:90,'Box','off')
ylabel('relative suitability for satisfyiing constraints')

print -dpdf -r350 suitabilityAllMeat.pdf
