load replaceAll2 C D X* Ipp

xbf = mean(Xbf,2);
xmt = mean(Xmt,2); clear X*

% g/cap-day
[xbf,ibf] = sort(xbf,'descend'); Ibf = Ipp(ibf);
[xmt,imt] = sort(xmt,'descend'); Imt = Ipp(imt);

% g protein/cap-day
protbf = xbf .* D(2,ibf)';
protmt = xmt .* D(2,imt)';

% land: m2/cap-day
landbf = xbf .* C(ibf,1);
landmt = xmt .* C(imt,1);

% gNR/cap-day
Nrbf = xbf .* C(ibf,2);
Nrmt = xmt .* C(imt,2);

% gCO2/cap-day
ghgbf = xbf .* C(ibf,3);
gngmt = xmt .* C(imt,3);

% Lwater/cap-day
watbf = xbf .* C(ibf,4);
watmt = xmt .* C(imt,4);  clear C D Ipp
