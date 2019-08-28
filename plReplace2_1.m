load replaceAll2 B* E* b* e* vrpp

% stats of replacement plant diets
mBbf = mean(Bbf,2);  sBbf = std(Bbf,0,2);
mBmt = mean(Bmt,2);  sBmt = std(Bmt,0,2);  clear B*

% find nutri constraints in terms of which beef or allMeat 
% normalized-differ markedly from their plant replacements
dbf = abs((mBbf - bbf)./bal); % beef   
dmt = abs((mBmt - bmt)./bal); % allMeat

% focus on significant diffs
i   = find( dbf>.25 & dmt>.25 );
h   = .7/length(i);

% normalized values, [ actual plantReplacement ] 
nv = 100*[   bmt(i)   bbf(i)  mBmt(i) mBbf(i) ]./(bal(i)*[1 1 1 1]);
ns = 100*[ zeros(length(i),2) sBmt(i) sBbf(i) ]./(bal(i)*[1 1 1 1]);
vr = vrpp(i);

[mx,i] = sort(max(nv,[],2));    clear mB* sB* vrpp
mx = .9*mx/max(mx);
nv = nv(i,:);
ns = ns(i,:);
vr = vr(i);
Np = length(i);

cl      = zeros(4,3);
cl(1,:) = [  0  0 .8 ]; % for all meat
cl(2,:) = [ .8  0  0 ]; % for beef
cl(3,:) = [  0  1  0 ]; % for all meat plant replacement
cl(4,:) = [  0 .7  0 ]; % for beef plant replacement

ha = 'HorizontalAlignment';
nr = 'Normalized';
lw = 'LineWidth';
ec = 'EdgeColor';
fs = 'FontSize';
un = 'Units';

clf

axes('Position',[.05 .05 .8 .8])
for i = 1:Np % loop on major nutri attributes
  y   = (i-1)*4.2;
  dh  = nv(i,:);
  dsh = ns(i,:);
  % the gray background
  if (i*.5)==floor(i*.5)
    x  = 170+(i>9)*105+(i==Np)*305;
    if i==Np; x = max(dh); end
    j  = fill([0 0 x x],y+[.1 4 4 .1],.9*[1 1 1]); hold on
    set(j,ec,.9*[1 1 1])
  end
  for j = 1:4
    fill([0 0 dh(j)*[1 1]],y+j-[.9 0 0 .9],cl(j,:),ec,cl(j,:)); hold on
    for k = 5:5:max(dh(j));  plot([k k],y+j-[0 .9],'w-',lw,.1);  end
  end  
  for j = 3:4
    if dh(j)<450 
      x = dh(j)+dsh(j)*.5*[-1 1];
      plot(x,y+j-.5*[1 1],'k-')
      plot(x(1)*[1 1],y+j-.5+.3*[-1 1],'k-')
      plot(x(2)*[1 1],y+j-.5+.3*[-1 1],'k-')
    else 
      x = dh(j)+dsh(j)*[-.5 0];
      plot(x,y+j-.5*[1 1],'k-')
      plot(x(1)*[1 1],y+j-.5+.3*[-1 1],'k-')
    end
  end
  if max(dh(:))<300;  text(max(dh(:)+.5*dsh(:))+2,y+3,vr(i)); 
  else;               text(418,y+1.4,vr(i),ha,'c'); end
end
set(gca,'YTick',[],'Box','off','XTick',0:25:max(nv(:)))
xlabel('percent of delivery by the truncated MAD')
axis([0 max(nv(:))+5 0 length(nv)*4.2+.5])
% legend
axes('Position',[.56 .53 .28 .13])
fill(   [0 0 8 8],9+[0 8 8 0],cl(2,:),ec,cl(2,:)); hold on
fill(   [0 0 8 8],  [0 8 8 0],cl(1,:),ec,cl(1,:))
fill(20+[0 0 8 8],9+[0 8 8 0],cl(4,:),ec,cl(4,:))
fill(20+[0 0 8 8],  [0 8 8 0],cl(3,:),ec,cl(3,:))
text( 9.2, 5.5,'all');
text( 9.2, 2.5,'meat')
text( 9.2,14.5,'beef')
text(29.2,14.5,'beef')
text(29.2,11.5,'replacement')
text(29.2, 5.5,'all meat')
text(29.2, 2.5,'replacement')
set(gca,'XTick',[],'YTick',[])
axis([-1 50 -1 18])

% stats of environmental costs of replacement plant diets
mEbf = mean(Ebf,2);  sEbf = std(Ebf,0,2);
mEmt = mean(Emt,2);  sEmt = std(Emt,0,2);

% next plots, environmental costs

t    = cell(4,1);
t{1} = 'm^2';
t{2} = 'g Nr';
t{3} = 'kg CO_{2eq}';
t{4} = 'L';
pn   = 'bcde';
x    = [ 1.3 2.3 1 2 ];

for i = 1:4 % loop on [ m2 gNr kgCO2e L ]/cap-day
  d = [ emt(i) ebf(i) mEmt(i) mEbf(i) ];
  s = [   0     0     sEmt(i) sEbf(i) ];
  p = round(100*d./eal(i));
  if i==3; d = d*1e-3; s = s*1e-3; end % if emissions, make it kg, not g
  j = max(d+s)*.05;
  axes('Position',[.36+(i-1)*.131 .075 .09 .35])
  if d(1)>d(3)
    fill(.7+.48*[-1 -1 1 1],[0 d(1)*[1 1] 0],cl(1,:),ec,cl(1,:)); hold on
    fill(.9+.48*[-1 -1 1 1],[0 d(3)*[1 1] 0],cl(3,:),ec,cl(3,:))
    plot([.9 .9],d(3)-s(3)*[1 0],'k-')	    
    plot([.9 .9],d(3)+s(3)*[1 0],'w-')	    
    plot(.9+.1*[-1 1],d(3)+s(3)*[1 1],'w-')
    plot(.9+.1*[-1 1],d(3)-s(3)*[1 1],'k-')
    text(.72,d(1)-j,[int2str(p(1)) '%'],ha,'c',fs,8,'Color','w')
    text(.92,d(3)-s(3)-j,[int2str(p(3)) '%'],ha,'c',fs,8)
  else
    fill(.7+.48*[-1 -1 1 1],[0 d(3)*[1 1] 0],cl(3,:),ec,cl(3,:)); hold on
    fill(.9+.48*[-1 -1 1 1],[0 d(1)*[1 1] 0],cl(1,:),ec,cl(1,:))    
    plot([.7 .7],d(3)+s(3)*[-1 1],'k-')	    
    plot(.7+.1*[-1 1],d(3)+s(3)*[1 1],'k-')
    plot(.7+.1*[-1 1],d(3)-s(3)*[1 1],'k-')
    text(.92,d(1)-j,[int2str(p(1)) '%'],ha,'c',fs,8,'Color','w')
    text(.72,d(3)-s(3)-j,[int2str(p(3)) '%'],ha,'c',fs,8)
  end    
  if d(2)>d(4)
    fill(2.0+.48*[-1 -1 1 1],[0 d(2)*[1 1] 0],cl(2,:),ec,cl(2,:))
    fill(2.2+.48*[-1 -1 1 1],[0 d(4)*[1 1] 0],cl(4,:),ec,cl(4,:))
    plot([2.2 2.2],d(4)+s(4)*[-1 1],'k-')	    
    plot(2.2+.1*[-1 1],d(4)+s(4)*[1 1],'k-')
    plot(2.2+.1*[-1 1],d(4)-s(4)*[1 1],'k-')
    text(2.02,d(2)-j,[int2str(p(2)) '%'],ha,'c',fs,8,'Color','w')
    text(2.22,d(4)-s(4)-j,[int2str(p(4)) '%'],ha,'c',fs,8)
  else
    fill(2.0+.48*[-1 -1 1 1],[0 d(4)*[1 1] 0],cl(4,:),ec,cl(4,:))
    fill(2.2+.48*[-1 -1 1 1],[0 d(2)*[1 1] 0],cl(2,:),ec,cl(2,:))
    plot([2 2],d(4)+s(4)*[-1 1],'k-')	    
    plot(2+.1*[-1 1],d(4)+s(4)*[1 1],'k-')
    plot(2+.1*[-1 1],d(4)-s(4)*[1 1],'k-')
    text(2.02,d(4)-s(4)-j,[int2str(p(4)) '%'],ha,'c',fs,8)
    text(2.22,d(2)-j,[int2str(p(2)) '%'],ha,'c',fs,8)
  end    
  if i==1
    text(1,1.14,'daily per capita environmental costs',un,nr)
  end
  text(.5,1.05,t{i},ha,'c',un,nr)
  set(gca,'XTick',[],'Box','off')
  text(-.08,1.04,['{\bf' pn(i) '}'],un,nr)
  axis([.08 2.68 0 max(d+s)*1.01])
end

print -dpdf -painters -r350 replace2_1.pdf %%%%%%%%%%%%%%%%%%%%%%%%


