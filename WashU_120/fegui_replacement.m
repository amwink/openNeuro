% add bias to the path
addpath(genpath('bias'));

% load fmris -- test both (warped) and (warped + smoothed)
funclist={'fd_cleanfuncs.txt'};

% test both 1 dynamic (so static) ECM and 100 dynamics
dynopts=[1 100]; 

% test both original (corr + 1) and relu (ReLU correlation)
coropts={'fast' 'relu'};

% set default options for fasteCM calls
opts.rankmap=0;
opts.normmap=0;
opts.degmap=0;
opts.maxiter=50;
opts.maskfile=[pwd filesep 'maskgroup.nii.gz'];
opts.atlasfile=0;
opts.wholemat=0;

for i=1:length(funclist)

  fid=fopen(funclist{i});
  funcfiles=textscan(fid,'%s');
  funcfiles=funcfiles{1};
  fclose(fid);

  for d=1:2
  
    opts.dynamics=dynopts(d);
    
    for c=1:2
      
      opts.correlate=coropts{c};
      
      for f=1:length(funcfiles)
	
	opts.inputfile=funcfiles{f};
	errorcode=fastECM(opts);
	
	if (~errorcode)
	  
	  dynstr=sprintf('dyn%03d',dynopts(d));
	  filespec=['fastECM' '_' coropts{c} '_' dynstr];
	  outfile=strrep(funcfiles{f},'.nii.gz','_fastECM.nii.gz');
	  newfile=strrep(outfile,'fastECM',filespec);
	  movefile(outfile,newfile);
	  
	end
	
      end % for f
      
    end % for c
    
  end % for d
  
end % for i
