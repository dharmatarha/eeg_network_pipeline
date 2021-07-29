% Surrigate Comparisons Wrapper script
%
%
%
%

freqs = {'delta', 'theta', 'alpha', 'beta', 'gamma'};
methods = {'plv', 'iplv', 'ciplv', 'ampCorr', 'orthAmpCorr'};
subjects = {'s02','s03','s04','s05','s06','s07','s08','s09',...
  's11','s12','s13','s14','s15','s16','s17','s18','s19','s20',...
  's21','s22','s23','s24','s25','s26','s27','s28'};
subjectsDelta = {'s02','s03','s04','s05','s06','s07','s08','s09',...
  's11','s12','s13','s14','s15','s16','s17','s18','s20',...
  's21','s22','s23','s24','s25','s26','s27'};


dirName = '/media/NAS502/adamb/hyperscan/newSurrEdgeEstimates/';


for f = 1:numel(freqs)
    freq = freqs{f};
    
    if strcmp(freq, 'delta')
        s = subjectsDelta;
    else
        s = subjects;
    end
    
    for m = 1:numel(methods)
        method = methods{m};
        
        cmpSurrRealConn_hyperscan4D(freq, dirName, s, method);
        
        cmpSurrRealConn_v2_hyperscan4D(freq, dirName, s, method);
        
    end
end

