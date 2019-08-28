% solve the diet optimization problem: min(x^T c) S.T. Ax >=< b 

i = '/Users/gidon/Dropbox (EshelMiloGroup)/PlantVmixedDiets/';

[D,T] = xlsread([i 'plVals2018.xlsx']);

Ic  = T(2:end,1); % all food items (plants AND animals)

% eliminate the end (Refs etc)
Ic  = Ic(1:length(D(:,1)));

% make variable names lowercase
vrp = lower(T(1,2:end)'); clear T

% make nutri info per g, not per 100 g
% each row in D is a food item
% each column in D is an attribute
D(:,1:54) = D(:,1:54)*.01; 

m = D(:,70); % dry to cooked mass ratio

% g/person-day of the iterms in the MAD
MADgf = D(:,67);               % kcal/person-day in MAD

i     = find(isnan(MADgf));    % if we have no info on MAD use of 
MADgf(i) = 1;                  % the item, assume it's vanishingly 
                               % small but nonzero, 1 kcal/person-day
MADgf = MADgf ./ D(:,1);       % kcal/cap-d * g/kcal = g/cap-day

% assemble the cost matrix: m2/g, gNr/g, gCO2e/g, L/g
C   =     D(:,61)*.001        ; % cropland needs, m2-yr/lossAdjustedGram
C   = [ C D(:,64)*453.592e-6 ]; % add Nr needs in gNr/lossAdjustedGram
C   = [ C D(:,73)            ]; % add GHG emissions in gCO2eq/lossAdjustedGram
C   = [ C D(:,76)*.001       ]; % add blue water needs in L/lossAdjustedGram

% from Alon: these are the values that do not need m correction:
% Choline:         column  27  ; row  49 
% vit k:           column  40  ; row  49
% phytosterol:     column  44  ; rows 21-24,49,58
% total flavonoid: column  51  ; rows 21,40,58
% soluble fiber:   column  52  ; rows 8,9,23,35,46,58,59,60,61,62,65
% omega3/6:        columns 53/4; rows 8,22,

% correct bulk dry items for cooking-added water except for above attributes
i = [ 1:26 28:39 41:43 45:50 55:length(vrp) ];  
D(:,i) = D(:,i) .* (m*ones(1,length(i)));

% correct choline & vit.K
j = [ 1:48 50:length(Ic) ];
D(j,[27 40]) = D(j,[27 40]) .* (m(j)*[1 1]);

% correct phytosterol 
j = [ 1:20 25:48 50:57 59:length(Ic) ];
D(j,44) = D(j,44) .* m(j);

% correct total falvonoid
j = [ 1:20 22:39 41:57 59:length(Ic) ];
D(j,51) = D(j,51) .* m(j);

% correct soluble fiber
j = [ 1:7 10:22 24:34 36:45 47:57 63:64 66:length(Ic) ];
D(j,52) = D(j,52) .* m(j);

% correct omega 3/6
j = [ 1:7 9:21 23:length(Ic) ];
D(j,[53 54]) = D(j,[53 54]) .* (m(j)*[1 1]);

% exclude columns (attributes) that we have either 
% already used (costs) or aren't going to
i   = [ 1:23 27:29 32:37 39:45 51:54 ];
D   = D(:,i);
vrp = vrp(i);

% set all remaining missing (exotic) values to zero
i = find(isnan(D)); D(i) = 0;

% make D be food items along rows, and nutri attributers along columns
% so now a row corresponds to a nutrtional attribute, 
% and a column corresponds to a food item
D = D'; 

% updated national resource usage

% cropland:  https://www.ers.usda.gov/ ...
%     data-products/major-land-uses/major-land-uses/#Cropland
nat = 392; % 2012 cropland, million acres

% reactive nitrogen in 2002 from all human activities in the
% US, from https://yosemite.epa.gov/sab/sabproduct.nsf/WebBOARD/INCFullReport
% ... /$File/Final%20INC%20Report_8_19_11(without%20signatures).pdf. 
% Table 1@ https://www.ers.usda.gov/data-products/fertilizer-use-and-price.aspx
% shows total primary N mass FERTILIZER use of 12.1e6 metric tons
nat = [ nat 28.5 ]; % Nr, Tg or million metric tons N

% 2017 GHG emissions, from https://www.epa.gov/ghgemissions/ ...
% inventory-us-greenhouse-gas-emissions-and-sinks-1990-2016
nat = [ nat 6457 ]; % GHG, million metric tons CO2e

% 2015 freshwater withdrawls from 
% https://pubs.usgs.gov/circ/1441/circ1441.pdf
i = 281e3;        % freshwater withdrawls, million gallons/day
i = i*365;        % freshwater withdrawls, million gallons/year
i = i*0.00378541; % freshwater withdrawls, million m3/year
nat = [ nat i ];  % freshwater withdrawls, million m3/year

% develop weights, the proportion of ALL animals in national totals
w(1) = 77.54+69.30; % total animal HQ cropland use ResourceFeedMain, Q12+Y12
w(2) = 13.4/2.2;    % total animal use of Nr, 1e9kg/yr (ResourceFeedMain AE15)
w(3) = 428.85763;   % mMtGHG/yr (PartioningPerCal H25)
w(4) = 47.4e3;      % total animal use of water, m3/yr (ResourceFeedMain AE17)
w    = w./nat;      % divide by total US use

% combined normalized costs, [m2, gNr, gCO2e, L]/g
Cc =  (w(1)*(C(:,1)'-mean(C(:,1)))./std(C(:,1))) ...
    + (w(2)*(C(:,2)'-mean(C(:,2)))./std(C(:,2))) ...  
    + (w(3)*(C(:,3)'-mean(C(:,3)))./std(C(:,3))) ...
    + (w(4)*(C(:,4)'-mean(C(:,4)))./std(C(:,4)));

D = [ D ; D(1,:)*0+1 ];
i = cell(length(vrp)+1,1);
i(1:length(vrp)) = vrp;
i{end} = 'total mass';
vrp = i;                

j = [ 69 73 72 76 ]; % indices of beef, pork, poultry, and grass-fed beef

bbf  = D(:,j([1 4]))*MADgf(j([1 4]));  % the RHS, beef
bmt  = D(:,j       )*MADgf(j       );  % the RHS, all meat
bal  = D            *MADgf          ;  % the RHS, everything

ebf  = C(j([1 4]),:)'*MADgf(j([1 4]));  % environ costs, beef
emt  = C(j       ,:)'*MADgf(j       );  % environ costs, all meat
eal  = C'            *MADgf          ;  % environ costs, everything

% exclude corn starch (21), and syrups (67-68)
i      = [ 1:20 22:66 69:length(Cc)  ]; 
C      = C(i,:);
Cc     = Cc(i)';
D      = D(:,i);
Ic     = Ic(i);
fMADgf = MADgf;
MADgf  = MADgf(i);
m      = m(i);

% retain only plants
i   = 1:65; 
Cf  = C;
Df  = D;
C   = C(i,:);
Cc  = Cc(i);
D   = D(:,i);
Ip  = Ic(i);
Np  = length(i); % number of plant items

% find garlic's index
for i = 1:Np; if sum(Ip{i}(1:4)=='GARL')==4; isGarlic = i; break; end; end

% pretty up plant food names (the last p is for pretty)
Ipp{ 1} = 'almonds';
Ipp{ 2} = 'apples';
Ipp{ 3} = 'apricots';
Ipp{ 4} = 'asparagus';
Ipp{ 5} = 'avocado';
Ipp{ 6} = 'banana';
Ipp{ 7} = 'barley';
Ipp{ 8} = 'kidney beans';
Ipp{ 9} = 'snap beans';
Ipp{10} = 'blueberries';
Ipp{11} = 'broccoli';
Ipp{12} = 'cabbage';
Ipp{13} = 'canola oil';
Ipp{14} = 'cantaloupe';
Ipp{15} = 'carrots';
Ipp{16} = 'cauliflower';
Ipp{17} = 'celery';
Ipp{18} = 'cherries';
Ipp{19} = 'collards';
Ipp{20} = 'corn flr';
Ipp{21} = 'corn grts';
Ipp{22} = 'swt corn';
Ipp{23} = 'cucumber';
Ipp{24} = 'garlic';
Ipp{25} = 'grapefruit';
Ipp{26} = 'grapes';
Ipp{27} = 'hazelnuts';
Ipp{28} = 'honeydew';
Ipp{29} = 'kiwi';
Ipp{30} = 'lemon';
Ipp{31} = 'lettuce';
Ipp{32} = 'macadamia';
Ipp{33} = 'oats';
Ipp{34} = 'olive oil';
Ipp{35} = 'onion';
Ipp{36} = 'oranges';
Ipp{37} = 'peaches';
Ipp{38} = 'peanuts';
Ipp{39} = 'pears';
Ipp{40} = 'green peas';
Ipp{41} = 'grn peppers';
Ipp{42} = 'pineapple';
Ipp{43} = 'pistachio';
Ipp{44} = 'potato';
Ipp{45} = 'pumpkin';
Ipp{46} = 'raspberries';
Ipp{47} = 'rice';
Ipp{48} = 'soybn oil';
Ipp{49} = 'spinach';    
Ipp{50} = 'smmr squash';
Ipp{51} = 'strawberries';
Ipp{52} = 'swt potato';
Ipp{53} = 'tomatoes';
Ipp{54} = 'walnuts';
Ipp{55} = 'watermelon';
Ipp{56} = 'wheat';
Ipp{57} = 'chickpeas';
Ipp{58} = 'lentils';
Ipp{59} = 'soybeans';
Ipp{60} = 'tofu';
Ipp{61} = 'buckwheat';   
Ipp{62} = 'sorghum';     
Ipp{63} = 'rye';      
Ipp{64} = 'spelt';
Ipp{65} = 'safflower oil';

vrpp{ 1} = 'kcal';		
vrpp{ 2} = 'protein';		
vrpp{ 3} = 'fat';		
vrpp{ 4} = 'carbs';		
vrpp{ 5} = 'total fiber';	
vrpp{ 6} = 'sugar';		
vrpp{ 7} = 'calcium';		
vrpp{ 8} = 'iron';		
vrpp{ 9} = 'magnesium';		
vrpp{10} = 'phosphorus';	
vrpp{11} = 'potassium';		
vrpp{12} = 'sodium';		
vrpp{13} = 'zinc';		
vrpp{14} = 'copper';		
vrpp{15} = 'manganese';		
vrpp{16} = 'selenium';		
vrpp{17} = 'vit C';		
vrpp{18} = 'thiamin';		
vrpp{19} = 'riboflavin';	
vrpp{20} = 'niacin';		
vrpp{21} = 'panto acid';	
vrpp{22} = 'vit B6';		
vrpp{23} = 'folate';		
vrpp{24} = 'choline';		
vrpp{25} = 'B12';		
vrpp{26} = 'vit A';	
vrpp{27} = '\alpha carot';	
vrpp{28} = '\beta carot';	
vrpp{29} = '\beta crypt';	
vrpp{30} = 'lycopene';		
vrpp{31} = 'lut+zea';		
vrpp{32} = 'vit E';		
vrpp{33} = 'vit D';		
vrpp{34} = 'vit K';		
vrpp{35} = 'FAsat';		
vrpp{36} = 'FAmono';		
vrpp{37} = 'FApoly';		
vrpp{38} = 'phytosterols';	
vrpp{39} = 'cholestrl';		
vrpp{40} = 'flavonoid';		
vrpp{41} = 'soluble fiber';	
vrpp{42} = '\Omega 3';		
vrpp{43} = '\Omega 6';		
vrpp{44} = 'mass';              

% get sineq, the vector whose elements are
% 1 for things we want low, -1 for ones we want high
sineq = ones(length(vrp),1); % 1 is lessThan, -1 for moreThan
i     = [ 2 5 7:11 13:34 36:38 40:43 ];  sineq(i) = -1;

%    <  1 energ_kcal         
%  2 >  2 protein_(g)	
%    <  3 lipid_tot_(g)	  
%    <  4 carbohydrt_(g)	  
%  5 >  5 fiber_td_(g)	
%    <  6 sugar_tot_(g)	  
%  7 >  7 calcium_(mg)	
%  8 >  8 iron_(mg)	  	
%  9 >  9 magnesium_(mg)	
% 10 > 10 phosphorus_(mg)	
% 11 > 11 potassium_(mg)	
%    < 12 sodium_(mg)	 	     
% 13 > 13 zinc_(mg)		
% 14 > 14 copper_(mg)		     
% 15 > 15 manganese_(mg)		     
% 16 > 16 selenium_(mcg)	
% 17 > 17 vit_c_(mg)		
% 18 > 18 thiamin_(mg)	
% 19 > 19 riboflavin_(mg)	
% 20 > 20 niacin_(mg)	
% 21 > 21 panto_acid_(mg)	
% 22 > 22 vit_b6_(mg)	
% 23 > 23 folate_tot_(mcg)	       
% 24 > 24 choline_tot_(mg)	       
% 25 > 25 vit_b12_(mcg)	       
% 26 > 26 vit_a_iu		       
% 27 > 27 alpha_carot_(mcg)	       
% 28 > 28 beta_carot_(mcg)	       
% 29 > 29 beta_crypt_(mcg)	       
% 30 > 30 lycopene_(mcg)	       
% 31 > 31 lut+zea_(mcg)	       
% 32 > 32 vit_e_(mg)		       
% 33 > 33 vit_d_iu		       
% 34 > 34 vit_k_(mcg)	       
%    < 35 fa_sat_(g)	 
% 36 > 36 fa_mono_(g)	       
% 37 > 37 fa_poly_(g)	       
% 38 > 38 phytosterols(mg) per 100 g    
%    < 39 cholestrl_(mg)	 
% 40 > 40 total flavonoid (mg/100g  
% 41 > 41 soluble fiber (g per 100g 
% 42 > 42 pfn311, total omega 3, gm 
% 43 > 43 pfn611, total omega 6, gm 
%    < 44 total mass                

save replaceAll2.mat bbf bmt bal w* nat C Cf D Df Ip* Ic vrp* *MADgf e* sineq

Nc = length(bbf);   clear *MADgf j m T % Nc is the # of constraints 

Ds = D .* (sineq*ones(1,Np));

Nmc = 500; % number of MC realizations
Nr  =  35; % number of plant items retained in each MC realization

% use the Bounds Method outlined by Paulraj S, Sumathi P. A
% Comparative Study of Redundant Constraints Identification Methods in
% Linear Programming Problems. Math Probl Eng. 2010;2010:1-16
% doi:10.1155/2010/723402 to determine constraint redundency

iub = find(sineq== 1); % upper bound constraints 
ilb = find(sineq==-1); % lower bound constraints (and protein)

ub =  ones(Np,1)*25;
% limit garlic's mass upper bound
ub(isGarlic) = 5;

Nst = zeros(Nc,2);
for i = 1:Nmc
  m   = sort(randperm(Np,35)); % retained plant items
  % lower bound constraints
  lb = max( 0 , randn(length(m),1)*.5 );
  for jj = 1:length(ilb)
    kk = ilb(jj);
    Nst(kk,1) = Nst(kk,1) + ((D(kk,m)*lb)>=bmt(kk)); % all-meat
    Nst(kk,2) = Nst(kk,2) + ((D(kk,m)*lb)>=bbf(kk)); % beef-only
  end                                       
  % upper bound constraints   
  for jj = 1:length(iub)                     
    kk = iub(jj);
    Nst(kk,1) = Nst(kk,1) + ((D(kk,m)*ub(m))<=bmt(kk)); % all-meat
    Nst(kk,2) = Nst(kk,2) + ((D(kk,m)*ub(m))<=bbf(kk)); % beef-only
  end
end
Nst  = 100*Nst/Nmc;            % times 100/500 to get percent
exc = find( sum(Nst,2)==200 ); % exclude due to redundancy

% %%%%%%

% exclusions: all redundent constraints above plus
% proein (equality; 2), carbs (4) & sugar (6), of which beef has none
% B12 (25) & vitD (33) of which plants have none
exc = unique([ exc' 2 4 6 25 33 ]);
jkc = 1:Nc;
jkc = jkc(~ismember(1:Nc,exc));

save -append replaceAll2.mat jkc Nst

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%      beef only replacement optimizations     %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = zeros(Np,1);
ub =  ones(Np,1)*25;
ub(isGarlic) = 5; % limit garlic's mass upper bound

% assign RHS, and relax some bounds
b      = bbf;
b(16)  = b(16)*0.5; % Se
b(24)  = b(24)*0.5; % Cholin
b(36)  = b(36)*0.1; % monoFA
b(44)  =       250; % total mass (instead of the actual 66)
bs     = b .* sineq;

X  = zeros(Np,Nmc);
B  = zeros(Nc,Nmc);
E  = zeros( 4,Nmc);      

Nh  = 0;    
Nal = 0;
while Nh < Nmc
  Nal = Nal+1;
  lbh = max( 0 , lb + randn(Np,1)*.5 );
  ubh =          ub + randn(Np,1)* 5  ;  
  ubh(isGarlic) = ub(isGarlic) + randn;
  m   = sort(randperm(Np,Nr)); % retained plant items
  lbh = lbh(m);
  ubh = ubh(m);  
  [x,~,e] = cplexlp(Cc(m),Ds(jkc,m),bs(jkc),D(2,m),b(2),lbh,ubh);
  if e==1 % see https://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.3/ ...
    % ilog.odms.cplexr.help/refmatlabcplex/html/cplexmiqcp-m.html
    % for cplexlp's exit flgs
    Nh      = Nh+1;
    X(m,Nh) = x;
    B(:,Nh) = D(:,m)*x;
    E(:,Nh) = C(m,:)'*x;
    disp(['beef only, ' int2str(Nh)])
  end
end
Bbf   = B;
Xbf   = X;
Ebf   = E;
Nalbf = Nal;

save -append replaceAll2.mat Bbf Ebf Xbf Nalbf

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%   all-Meat replacement optimizations   %%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lb = zeros(Np,1);
ub =  ones(Np,1)*25;
ub(isGarlic) = 5; % limit garlic's mass upper bound

% assign RHS, and relax some bounds
b     = bmt;
b(16) = b(16)*0.5; % Se
b(20) = b(20)*0.5; % Niacin (B3)
b(24) = b(24)*0.5; % Cholin
b(36) = b(36)*0.5; % monoFA
b(44) =       300; % total mass (instead of the actual 256)

bs    = b .* sineq;

X  = zeros(Np,Nmc);
B  = zeros(Nc,Nmc);
E  = zeros( 4,Nmc);      

Nh  = 0;    
Nal = 0;
while Nh < Nmc
  Nal = Nal+1;
  lbh = max( 0 , lb + randn(Np,1)*.5 );
  ubh =          ub + randn(Np,1)* 5  ;  
  ubh(isGarlic) = ub(isGarlic) + randn;
  m   = sort(randperm(Np,Nr)); % retained plant items
  lbh = lbh(m);
  ubh = ubh(m);  
  [x,~,e] = cplexlp(Cc(m),Ds(jkc,m),bs(jkc),D(2,m),b(2),lbh,ubh);
  if e==1 % see https://www.ibm.com/support/knowledgecenter/SSSA5P_12.6.3/ ...
    % ilog.odms.cplexr.help/refmatlabcplex/html/cplexmiqcp-m.html
    % for cplexlp's exit flgs
    Nh      = Nh+1;
    X(m,Nh) = x;
    B(:,Nh) = D(:,m)*x;
    E(:,Nh) = C(m,:)'*x;
    disp(['all meat, ' int2str(Nh)])
  end
end
Bmt   = B;
Xmt   = X;
Emt   = E;
Nalmt = Nal;

save -append replaceAll2.mat Bmt Emt Xmt Nalmt
