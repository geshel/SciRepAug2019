load replaceAll2 B* b* vrpp sineq Nst

Nc = length(Bbf(:,1));

% diet composition differences: stats over replacement diets
Dbf = Bbf - bbf;
Dmt = Bmt - bmt;   

mbf = mean(Bbf,2);
mmt = mean(Bmt,2);

mdbf = mean(Dbf,2);   sdbf = std(Dbf,0,2);
mdmt = mean(Dmt,2);   sdmt = std(Dmt,0,2);

disp('--------------------------------------------------------------')
disp('          cricical inequality constraints, beef')
disp('--------------------------------------------------------------')
disp('        |          |       |      | normalized |        |')
disp('        |   mean   |  RHS  | mean |    mean    |        |')
disp('item    | solution | bound | diff |    diff    | 100s/m | Nst')
disp('--------------------------------------------------------------')

v = 100*sdbf./mbf;
n = 100*mdbf./mbf;
j = find( (n.*(-sineq))<5 & Nst(:,2)<5 ); 
for k = 1:length(j)
  i  = j(k); 
  tt = [ vrpp{i} char(' '*ones(1,13-length(vrpp{i}))) ];
  tt = [ tt sprintf('%6.1f',mbf(i)) ' ' ];
  if sineq(i)==1; sg = '<'; else; sg = '>'; end
  tt = [ tt sg sprintf('%5.1f', bbf(i)) ' ' ];
  tt = [ tt    sprintf('%6.1f',mdbf(i)) ' ' ];
  tt = [ tt    sprintf('%9.1f',   n(i)) ' ' ];
  tt = [ tt    sprintf('%9.1f',   v(i)) ' ' ];
  tt = [ tt    sprintf('%7.1f', Nst(i,2))   ];
  disp(tt) 
end

disp(' ')
disp('------------------------------------------------------------')
disp('         cricical inequality constraints, all-meat')
disp('------------------------------------------------------------')

v = 100*sdmt./mmt;
n = 100*mdmt./mmt;
j = find( (n.*(-sineq))<5 & Nst(:,1)<5 );
for k = 1:length(j)
  i  = j(k); 
  tt = [ vrpp{i} char(' '*ones(1,13-length(vrpp{i}))) ];
  tt = [ tt sprintf('%6.1f',mmt(i)) ' ' ];
  if sineq(i)==1; sg = '<'; else; sg = '>'; end
  tt = [ tt sg sprintf('%5.1f', bmt(i)) ' ' ];
  tt = [ tt    sprintf('%6.1f',mdmt(i)) ' ' ];
  tt = [ tt    sprintf('%9.1f',   n(i)) ' ' ];
  tt = [ tt    sprintf('%9.1f',   v(i)) ' ' ];
  tt = [ tt    sprintf('%7.1f', Nst(i,1))   ];
  disp(tt) 
end

disp('------------------------------------------------------------')
