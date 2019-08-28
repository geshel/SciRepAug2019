
load replaceAll2

Np  = length(Xmt(:,1)); % plant item num
Nmc = length(Xmt(1,:)); % MC realization num
No  = length(  C(1,:)); % environ cost num 
Nc  = length(  D(:,1)); % number of constraints (active or not)

% plant item mass stats in the animal replacement diets
Xmmt = mean(Xmt,2);   Xsmt = std(Xmt,0,2);
Xmbf = mean(Xbf,2);   Xsbf = std(Xbf,0,2);

% diet composition stats in the animal replacement diets
Bmmt = mean(Bmt,2);   Bsmt = std(Bmt,0,2);
Bmbf = mean(Bbf,2);   Bsbf = std(Bbf,0,2);

d     = bmt - bbf; % diff in content per g
d     = 100*d./bal;     % normalize to percent of MAD delivery  
[~,i] = sort(abs(d),'descend');

disp('------------------------------------------------------------')
disp(' cricical active inequality constraints, beef')
disp('--------------------------------------------------------------')
disp('          |          |       |      | normalized |        ')
disp('          |   mean   |  RHS  | mean |    mean    |        ')
disp('item      | solution | bound | diff |    diff    | 100s/m ')
disp('--------------------------------------------------------------')

cbf = []; % for critical constraint indices

d = 100*Bsbf./Bmbf;
m = Bmbf-bbf;
n = 100*m./Bmbf;
[~,j] = sort(abs(n));
for i = 1:Nc
  k = j(i);
  if abs(n(k))<15 & d(k)<15
    tt = [ vrpp{k} char(' '*ones(1,15-length(vrpp{k}))) ];
    tt = [ tt sprintf('%6.1f',Bmbf(k)) ' ' ];
    if sineq(k)==1; sg = '<'; else; sg = '>'; end
    tt = [ tt sg sprintf('%5.1f',bbf(k)) ' ' ];
    tt = [ tt    sprintf('%6.1f',  m(k)) ' ' ];
    tt = [ tt    sprintf('%9.1f',  n(k)) ' ' ];
    tt = [ tt    sprintf('%9.1f',  d(k))     ];
    disp(tt) 
    cbf = [ cbf k ];
  end
end
cbf = unique(cbf);

disp(' ')
disp('------------------------------------------------------------')
disp(' cricical active inequality constraints, all-meat')
disp('------------------------------------------------------------')

cmt = []; % for critical constraint indices

d = 100*Bsmt./Bmmt;
m = Bmmt-bmt;
n = 100*m./Bmmt;
[~,j] = sort(abs(n));
for i = 1:Nc
  k = j(i);
  if abs(n(k))<15 & d(k)<15
    tt = [ vrpp{k} char(' '*ones(1,15-length(vrpp{k}))) ];
    tt = [ tt sprintf('%6.1f',Bmmt(k)) ' ' ];
    if sineq(k)==1; sg = '<'; else; sg = '>'; end
    tt = [ tt sg sprintf('%5.1f',bmt(k)) ' ' ];
    tt = [ tt    sprintf('%6.1f',  m(k)) ' ' ];
    tt = [ tt    sprintf('%9.1f',  n(k)) ' ' ];
    tt = [ tt    sprintf('%9.1f',  d(k))     ];
    disp(tt) 
    cmt = [ cmt k ];
  end
end
cmt = unique(cmt);

disp('------------------------------------------------------------')
disp(' ')
disp(' ')

Dmbf = [];
Dmmt = [];
for i = 1:length(cbf)-1 % exclude total mass, bcz its entries are all 1g/g
  d  = D(cbf(i),:)';
  d  = (d-mean(d))/std(d);
  % if we want the least, d<0 is desireable, so
  if sineq(cbf(i))==1; d = -d; end 
  Dmbf = [ Dmbf d ];  
end
dmbf = sum(Dmbf,2);
[dmbf,i] = sort(dmbf,'descend'); Imbf = Ipp(i(1:10));

Dmmt = [];
for i = 1:length(cmt)-1 % exclude total mass, bcz its entries are all 1g/g
  d  = D(cmt(i),:)';
  d  = (d-mean(d))/std(d);
  % if we want the least, d<0 is desireable, so
  if sineq(cmt(i))==1; d = -d; end 
  Dmmt = [ Dmmt d ];  
end
dmmt = sum(Dmmt,2);
[dmmt,i] = sort(dmmt,'descend'); Immt = Ipp(i(1:10));

Nscb = Nmc/Nalbf; % fraction of feasible among total MC trials, beef-only 
Nscb = round(Nscb*1e6);  % successful MC cases / million tries, beef-only 
Nsca = Nmc/Nalmt; % fraction of feasible among total MC trials, all-meat 
Nsca = round(Nsca*1e6);  % successful MC cases / million tries, all-meat 

disp(' ')
disp('------------------------------------------------------------')
i = 'successful MC realizations per 1e6 trials,  all meat = ';
disp([ i num2str(Nsca) ])
i = 'successful MC realizations per 1e6 trials, beef only = ';
disp([ i num2str(Nscb) ])
disp('------------------------------------------------------------')
disp(' ')

clear all %%%% start fresh %%%%%%%%%%

load replaceAll2

Np  = length(Xmt(:,1)); % plant item num
Nmc = length(Xmt(1,:)); % MC realization num
No  = length(  C(1,:)); % environ cost num 
Nc  = length(  D(:,1)); % number of constraints (active or not)

% plant item mass stats in the animal replacement diets
Xmmt = mean(Xmt,2);   Xsmt = std(Xmt,0,2);
Xmbf = mean(Xbf,2);   Xsbf = std(Xbf,0,2);

% diet composition stats in the animal replacement diets
Bmmt = mean(Bmt,2);   Bsmt = std(Bmt,0,2);
Bmbf = mean(Bbf,2);   Bsbf = std(Bbf,0,2);

% start Erratum-writing tests %%%%%%%%%%%%%%%%

% item specific protein contributions to total protein in the total diet
i = D(2,:)' * ones(1,Nmc);
PRmt  = Xmt.*i;
PRbf  = Xbf.*i;
spmt  = ones(Np,1)*sum(PRmt,1);
spbf  = ones(Np,1)*sum(PRbf,1);
nPRmt = 100*mean((PRmt./spmt),2);
nPRbf = 100*mean((PRbf./spbf),2);

% item specific absolute environmental costs
Eibf = zeros(Np,No,Nmc);     Eimt = Eibf;    j = ones(1,Nmc);          
i = C(:,1)*j; Eibf(:,1,:) = Xbf.*i; Eimt(:,1,:) = Xmt.*i; 
i = C(:,2)*j; Eibf(:,2,:) = Xbf.*i; Eimt(:,2,:) = Xmt.*i; 
i = C(:,3)*j; Eibf(:,3,:) = Xbf.*i; Eimt(:,3,:) = Xmt.*i; 
i = C(:,4)*j; Eibf(:,4,:) = Xbf.*i; Eimt(:,4,:) = Xmt.*i; 
mEibf = mean(Eibf,3);    sEibf = std(Eibf,[],3);  
mEimt = mean(Eimt,3);    sEimt = std(Eimt,[],3);    

% normalize item specific environ costs to % of respective costs' totals
tEbf = squeeze(sum(Eibf,1));   pEbf = zeros(size(Eibf));
tEmt = squeeze(sum(Eimt,1));   pEmt = zeros(size(Eimt));   j = ones(Np,1);
for i = 1:No
  pEbf(:,i,:) = 100*squeeze(Eibf(:,i,:))./(j*tEbf(i,:));   
  pEmt(:,i,:) = 100*squeeze(Eimt(:,i,:))./(j*tEmt(i,:));   
end
mpEbf = mean(pEbf,3);    spEbf = std(pEbf,[],3);  
mpEmt = mean(pEmt,3);    spEmt = std(pEmt,[],3);    

% which items dominate over all environ costs?
j = 1:7;
[~,i1] = sort(mpEbf(:,1),'descend');
[~,i2] = sort(mpEbf(:,2),'descend');
[~,i3] = sort(mpEbf(:,3),'descend');
[~,i4] = sort(mpEbf(:,4),'descend'); 
ib = unique([ i1(j) ; i2(j) ; i3(j) ; i4(j) ]); 
[~,i1] = sort(mpEmt(:,1),'descend');
[~,i2] = sort(mpEmt(:,2),'descend');
[~,i3] = sort(mpEmt(:,3),'descend');
[~,i4] = sort(mpEmt(:,4),'descend'); 
im = unique([ i1(j) ; i2(j) ; i3(j) ; i4(j) ]); 
j  = unique([im;ib]);  clear i*
xmmmt = Xmmt(j);
xmmmt = Xmmt(j);
xmmbf = Xmbf(j);
xmmbf = Xmbf(j);
protmt = nPRmt(j);
protbf = nPRbf(j);
envmt  = mpEmt(j,:);
envbf  = mpEbf(j,:);
itnms  = Ipp(j); 

% end of Erratum-writing tests %%%%%%%%%%%%%%%%

% item specific environmental contributions
atEmt = zeros(   No,Nmc); % absolute total, meat
atEbf = zeros(   No,Nmc); % absolute total, beef
MEmt  = zeros(   No,Nmc); % percent of the replaced diet, meat
MEbf  = zeros(   No,Nmc); % percent of the replaced diet, beef
Emt   = zeros(Np,No,Nmc); % item specific percent contribution to total, meat
Ebf   = zeros(Np,No,Nmc); % item specific percent contribution to total, beef
for i = 1:Nmc
  E1 = C .* (Xmt(:,i)*[1 1 1 1]);  s1 = sum(E1);
  E2 = C .* (Xbf(:,i)*[1 1 1 1]);  s2 = sum(E2);  
  atEmt(:,i) = s1';
  atEbf(:,i) = s2';
  MEmt(:,i) = 100*s1'./emt;
  MEbf(:,i) = 100*s2'./ebf;
  % normalize to percent of each cost sum
  Emt(:,:,i) = 100*E1./(ones(Np,1)*s1);
  Ebf(:,:,i) = 100*E2./(ones(Np,1)*s2);
end

matEmt = mean(atEmt,2); satEmt = std(atEmt,0,2); 
matEbf = mean(atEbf,2); satEbf = std(atEbf,0,2);  clear atE* E1 E2

% absolute annual national envir costs, w units the same as nat:
% cropland, million acres
% Nr, Tg or million metric tons N
% GHG, million metric tons CO2e
% fresh water withdrawl, million m3
% (recal that C units are m2/g, gNr/g, gCO2e/g, L/g)
toNat = 365*327.* [ (1/4046.86) 1e-6 1e-6 1e-3 ];
matEmt = matEmt'.*toNat;
matEbf = matEbf'.*toNat;
satEmt = satEmt'.*toNat;
satEbf = satEbf'.*toNat;

mEmt   = mean( MEmt,2);   sEmt = std( MEmt,0,2);
mEbf   = mean( MEbf,2);   sEbf = std( MEbf,0,2);
Emmt   = mean(  Emt,3);   Esmt = std(  Emt,0,3);
Embf   = mean(  Ebf,3);   Esbf = std(  Ebf,0,3);  clear Emt Ebf ME* s1 s2

prmt = D(2,:) .* Xmmt'; % protein in the mean diet, all meat
prbf = D(2,:) .* Xmbf'; % protein in the mean diet, beef

prmt = 100*prmt/bmt(2); % percent of total protein, all meant
prbf = 100*prbf/bbf(2); % percent of total protein, beef 

% find out how many food items are in the top np items of all costs
np  = 7; 
Ncl = [];
for ic = 1:No
  [~,i] = sort(Embf(:,ic),'descend');  Ncl = [ Ncl i(1:np)' ];
  [~,i] = sort(Emmt(:,ic),'descend');  Ncl = [ Ncl i(1:np)' ];
end 
icl  = unique(Ncl);
Xmbf = Xmbf(icl,:);  Xmmt = Xmmt(icl,:);
Xsbf = Xsbf(icl,:);  Xsmt = Xsmt(icl,:);
Embf = Embf(icl,:);  Emmt = Emmt(icl,:);
Esbf = Esbf(icl,:);  Esmt = Esmt(icl,:);
prbf = prbf(icl  );  prmt = prmt(icl  );

% color map; start with too many colors
cl2 = [ 166 206 227 ; 106  61 154 ;  77 146  33 ; 222 119 174 ; ...
         39 100  25 ; 197  27 125 ; 253 224 239 ; 227  26  28 ; ...
	127 188  65 ;  51 160  44 ; 255 255 153 ; 177  89  40 ; ...
	 31 120 180 ; 241 182 218 ; 178 223 138 ; 202 178 214 ; ...
	230 245 208 ; 253 191 111 ; 184 225 134 ; 255 127   0 ; ...
	142   1  82 ; 251 154 153 ];
cl2 = cl2(1:length(icl),:)/255;

v      = 'mtbf';
clr    = [ 234 102 0 ; 193 0 76 ; 2 144 50 ; 92 92 215 ] / 255;
	 
rpl    = cell(2,1);
rpl{1} = 'all meat'; 
rpl{2} = 'beef only';

cst    = cell(4,1);
cst{1} = 'land';
cst{2} = 'N fertilier';
cst{3} = 'GHG \uparrow';
cst{4} = 'water';

pn = 'bcdefghijk';

ha = 'HorizontalAlignment';
va = 'VerticalAlignment';
nr = 'normalized';
ec = 'EdgeColor';
lw = 'LineWidth';
fs = 'FontSize';
un = 'Units';

clf

[~,ip] = sort(prmt+prbf,'descend');

% the food item color legend
axes('Position',[.01 .15 .085 .83])
set(gca,'Visible','off')
for i = ip
  j = ip(i);
  fill([1 1 8 8],(i-1)*10+[0 8 8 0],cl2(j,:),ec,cl2(j,:)); hold on
  Ih = Ipp{icl(j)};
  k = find(Ih(:)==' ');
  if isempty(k)
    if length(Ih(:))>8
      text(4.5,(i-1)*10+6.7,[Ih(1:6) '-'],ha,'c',fs,9)
      text(4.5,(i-1)*10+3.3,['-' Ih(7:end)],ha,'c',fs,9)
    else
      text(4.5,(i-1)*10+6.0,Ih,ha,'c',fs,9)
    end
  else
    text(4.5,(i-1)*10+6.7,Ih(1:k-1),ha,'c',fs,9)
    text(4.5,(i-1)*10+3.3,Ih(k+1:end),ha,'c',fs,9)
  end
end
set(gca,'Visible','off')

% the protein bottom axis
axes('Position',[.096 .14 .1164 .04])
plot([0 30],[1 1],'k-'); hold on
for i = 0:10:30
  plot([i i],[1 3],'k-')
  text(i+.2,-4.5,int2str(i),ha,'c')
end
text(-8,-13,'% of total protein')
axis([0 max([prbf prmt]) 0 9])
set(gca,'Visible','off','XLim',[0 30])

% the actual protein bars
axes('Position',[.096 .15 .1164 .83])
for i = 1:length(ip)
  j = ip(i);
  fill([0 0 prmt(j)*[1 1]],(i-1)*10+[4.2 8 8 4.2],cl2(j,:),ec,cl2(j,:)); hold on
  fill([0 0 prbf(j)*[1 1]],(i-1)*10+[0 3.8 3.8 0],cl2(j,:),ec,cl2(j,:))
  for k = 5:5:prmt(j)
    plot([k k],(i-1)*10+[4.2 8],'-','Color',.9*[1 1 1],lw,.2)
  end
  for k = 5:5:prbf(j)
    plot([k k],(i-1)*10+[0 3.8],'-','Color',.9*[1 1 1],lw,.2)
  end
end
set(gca,'Visible','off')
yr = get(gca,'YLim');
axis([0 30 yr])
text(-1,yr(2)+2,'{\bfa}')
text(6,yr(2)+2,'{\bfprotein}',fs,12)
text(2,yr(2)-13.1,'replacing meat',fs,8)
text(3,yr(2)-17.5,'replacing beef',fs,8)

for iv = 1:2 % loop on meat/beef
  % find top np env cost plant items in the animal replacement diet
  eval(['Xmh = Xm' v((iv-1)*2+[1 2]) ';'])
  eval(['Emh = Em' v((iv-1)*2+[1 2]) ';'])
  eval(['Xsh = Xs' v((iv-1)*2+[1 2]) ';'])
  eval(['Esh = Es' v((iv-1)*2+[1 2]) ';'])
  eval(['mEh = mE' v((iv-1)*2+[1 2]) ';'])
  eval(['sEh = sE' v((iv-1)*2+[1 2]) ';'])
  eval(['prh = pr' v((iv-1)*2+[1 2]) ';'])
  for ic = 1:4   % loop on cost: land, Nr, GHG, water
    [~,i] = sort(Emh(:,ic),'descend');
    i     = i(1:np);
    ph    = prh(i);
    clh   = cl2(i,:);
    Ih    = Ipp(icl(i));
    xmh   = Xmh(i);   
    xsh   = Xsh(i);
    emh   = Emh(i,ic);  
    esh   = Esh(i,ic);
    % now arrange these items, which contribute the most to the env
    % cost in question, in descending mass order
    [~,i] = sort(xmh,'descend');
    ph    = ph(i);
    clh   = clh(i,:);
    Ih    = Ih(i);
    xmh   = xmh(i);   
    xsh   = xsh(i);
    emh   = [ 0 ; cumsum(emh(i)) ];  
    esh   = [ 0 ;        esh(i)  ];     
    axes('Position',[.27+(ic-1)*.171 .98-iv*.43 .168 .42])
    hold on
    for i = 1:np
      fill([0 0 xmh(i)*[1 1]],emh([i i+1 i+1 i]),clh(i,:),ec,'w')
      if ic>3
        text(1.07,.5,[ 'replacing ' rpl{iv} ],un,nr,ha,'c','Rotation',90)
      end
      if iv<2; set(gca,'XTickLabel',[]); end
      if emh(i+1)>8
	x1 = 1;          y1 = emh(i)+2;
	x2 = x1+xsh(i);  y2 = emh(i)+2+esh(i+1);
        plot([x1 x2],[y1 y1],'w-')
        plot([x1 x1],[y1 y2],'w-')
      end
    end;  clear x1 x2 y1 y2
    for i = 5:5:40
      plot([i i],[-1 .6],'k-')
    end
    for i = 10:10:90
      plot([-1 .4],[i i],'k-')
    end
    if iv==1
      text(.5,1.06,['{\bf' cst{ic} '}'],un,nr,ha,'c',fs,12)
    end
    i = [ int2str(round(mEh(ic))) '\pm' ];
    i = [ i int2str(round(sEh(ic))) '%'];
    text(.52,.99,i,un,nr,ha,'c',fs,9)
    text(.98,.98,[ '{\bf' pn((iv-1)*4+ic) '}'],un,nr,ha,'r')
    axis([-.8 29 -1 88])
    if iv==2 & ic==2
      i = '                              ';
      xlabel([ i ' plant item mass, g (person\times d)^{-1}'])
    end
    i = '                             ';
    if ic>1; set(gca,'YTickLabel',[]); else; 
    if iv>1; ylabel([i i 'cumulative percent of resource use']); end; end
  end
  if iv==1 & ic==4
  end
end

print -dpdf -r350 replace2_3.pdf 
