function alignc = balign_peps(seqs, lea, lepmin, lepmax, weights) % lepmin lepmax indicate what lengths to include in the msa on which HMM is learnt
yes_red = 0;

yes_weight = 0;
if nargin > 4 & weights
 weights = importdata('weights.txt');
 Weights = weights';
 yes_weight = 1;
end

Xdseq=seqs;
M=length(seqs);

if size(Xdseq) == [1 M],
 Xdseq = reshape(seqs, [M 1]);
end

l=zeros(M,1);
lold=0;
  for i=1:M, 
     l(i)=size(Xdseq{i,1},2);   % curl brackets matrix of seqs of different length; for this object need to use cell
  end
lmax=max(l); 
lmin=min(l);

if lmin==lmax,
	af = cell2mat(Xdseq);
	N = length(af(1,:));
	alignc = zeros(M,N);
	for i=1:M
	    for j=1:N
		%same alphabet as matlab seqprofile used here in letter2number
	       alignc(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
	    end
	end
return
end

if lepmax > lmax
lepmax = lmax;
end
if lepmin < lmin
lepmin = lmin;
end

indexl=cell(lmax,1);
dist=zeros(lmax,1);
for ll=1:lmax
 indexl{ll}=find(l==ll); % indices of seqs. in each length group
 dist(ll)=length(indexl{ll}); % population of each length group
end                
al=cell(lmax,1);
p=cell(lmax,1);
W=cell(lmax,1);
logo=cell(lmax,1);
for ll=lmin:lmax 
al{ll} = cell2mat(Xdseq(indexl{ll},1)); % converts a cell array into an ordinary array
if yes_weight
W{ll} = Weights(indexl{ll});
p{ll} = seqprofile_rew(al{ll},21,W{ll});
else
p{ll}  = seqprofile(Xdseq(indexl{ll},1),'gaps','all','counts',true); % seq profiles (i.e. aa count in each position) for each length
%p{ll}  = seqprofile(Xdseq(indexl{ll},1),'gaps','all'); % seq profiles (i.e. aa count in each position) for each length
end
end

nl=0; 
lo=lepmin; 
a_old=al{lepmin,1};
a_new2 = a_old;
l_old=dist(lepmin); 
p_old= p{lepmin};
if yes_weight
w_old = W{lepmin};
end
if lepmin < lepmax
if lepmin + 1 < lepmax % lep decides which lengths should be considered in the alignment
for ll=lepmin + 1 : lepmax
    nl=nl+1;
    if (dist(ll)>0)
	    [pp,h1,h2] = profalign(p_old, p{ll}, 'ScoringMatrix', 'BLOSUM62');
	    a_new = repmat('-',l_old+dist(ll),ll);
	    a_new(1:l_old,h1) = a_old;
	    a_new(l_old+1:l_old+dist(ll),h2)=al{ll,1};
 	    if yes_weight
            w_new = [w_old, W{ll}];
	    end
    else
     a_new=a_old; 
     if yes_weight
     w_new = w_old;  
     end
    end
    if ll == lea
      a_new2 = a_new;
    end
    a_old = a_new;
    l_old=size(a_new,1); 
    if yes_weight
     w_old = w_new;
    p_old = seqprofile_rew(a_old,21,w_old);
    else
     p_old=seqprofile(a_old,'gaps','all','counts',true); % contains all previous lengths
    end
end
elseif lepmin==lea
    ll=lepmin;
    nl=nl+1;
    if(dist(ll)>0)
    [pp,h1,h2]=profalign(p_old, p{ll}, 'ScoringMatrix', 'BLOSUM62');
    a_new = repmat('-', l_old+dist(ll),ll);
    a_new(1:l_old,h1) = a_old;
    a_new(l_old+1:l_old+dist(ll),h2)=al{ll,1};
    if yes_weight
            w_new = [w_old, W{ll}];
	    end
    else
     a_new=a_old;  
     if yes_weight
     w_new = w_old;  
     end 
    end
    if ll == lea
      a_new2 = a_new;
    end
    a_old = a_new;
    l_old=size(a_new,1); 
     if yes_weight
     w_old = w_new;
    p_old = seqprofile_rew(a_old,21,w_old);
    else
     p_old=seqprofile(a_old,'gaps','all','counts',true); % contains all previous lengths
    end
elseif lepmax == lea
  ll=lepmax;
    nl=nl+1;
    if(dist(ll)>0)
    [pp,h1,h2]=profalign(p_old, p{ll}, 'ScoringMatrix', 'BLOSUM62');
    a_new = repmat('-', l_old+dist(ll),ll);
    a_new(1:l_old,h1) = a_old;
    a_new(l_old+1:l_old+dist(ll),h2)=al{ll,1};
    if yes_weight
            w_new = [w_old, W{ll}];
	    end
    else
     a_new=a_old;  
     if yes_weight
     w_new = w_old;  
     end 
    end
    if ll == lea
      a_new2 = a_new;
    end
    a_old = a_new;
    l_old=size(a_new,1); 
     if yes_weight
     w_old = w_new;
    p_old = seqprofile_rew(a_old,21,w_old);
    else
     p_old=seqprofile(a_old,'gaps','all','counts',true); % contains all previous lengths
    end
end
else
a_new = a_old;
a_new2 = a_old
end

Model = hmmprofstruct(lea);
hmmodel = hmmprofestimate(Model, a_new);

if lea < lmax
if lea + 1 < lmax
  for ll = lea + 1 : lmax
    if(dist(ll)>0)
    newseq = al{ll,1}; 
    for ns=1:dist(ll)
        [scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
        [zz,indal]=hmmprofmerge(aligned_seqs);
        ali_seq = aligned_seqs(indal);
	if length(newseq(ns,:)) ~= 9 & length(strfind(ali_seq,'-')) > 1 & yes_red
		if length(newseq(ns,:)) == 8
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)-1
		 seqtt = [seqt(1:t), '-', seqt(t+1:end)];
		 allist= [allist; seqtt];
		 scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 10
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 11
		seqt = newseq(ns,:);
		scorelist = [];
		allist=[];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		for s=1:length(seqtt)
			seqt3 = seqtt;
			seqt3(s) = [];
			allist = [allist; seqt3];
			scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end
        a_new2=[a_new2;ali_seq];
    end

    end
  end
 else
    ll = lmax;
    if(dist(ll)>0)
    newseq=al{ll,1}; 
    for ns=1:dist(ll)
        [scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
        [zz,indal]=hmmprofmerge(aligned_seqs);
        ali_seq=aligned_seqs(indal);
	if length(newseq(ns,:)) ~= 9 & length(strfind(ali_seq,'-')) > 1 & yes_red
		if length(newseq(ns,:)) == 8
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)-1
		 seqtt = [seqt(1:t), '-', seqt(t+1:end)];
		 allist= [allist; seqtt];
		 scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 10
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 11
		seqt = newseq(ns,:);
		scorelist = [];
		allist=[];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		for s=1:length(seqtt)
			seqt3 = seqtt;
			seqt3(s) = [];
			allist = [allist; seqt3];
			scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end
        a_new2=[a_new2;ali_seq];
    end
    end
  end
end
if lmin < lepmin
if lmin < lepmin - 1
  for ll = lmin:lepmin-1
   if(dist(ll)>0)
    newseq = al{ll,1}; 
    for ns=1:dist(ll)
        [scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
        [zz,indal]=hmmprofmerge(aligned_seqs);
        ali_seq=aligned_seqs(indal);

	if length(newseq(ns,:)) ~= 9 & length(strfind(ali_seq,'-')) > 1 & yes_red
		if length(newseq(ns,:)) == 8
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)-1
		 seqtt = [seqt(1:t), '-', seqt(t+1:end)];
		 allist= [allist; seqtt];
		 scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 10
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 11
		seqt = newseq(ns,:);
		scorelist = [];
		allist=[];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		for s=1:length(seqtt)
			seqt3 = seqtt;
			seqt3(s) = [];
			allist = [allist; seqt3];
			scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end

        a_new2=[a_new2;ali_seq];
     end
    end
  end
 else
    ll=lmin;
    if(dist(ll)>0)
    newseq=al{ll,1}; 
    for ns=1:dist(ll)
        [scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
        [zz,indal]=hmmprofmerge(aligned_seqs);
        ali_seq=aligned_seqs(indal);

	if length(newseq(ns,:)) ~= 9 & length(strfind(ali_seq,'-')) > 1 & yes_red
		if length(newseq(ns,:)) == 8
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)-1
		 seqtt = [seqt(1:t), '-', seqt(t+1:end)];
		 allist= [allist; seqtt];
		 scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 10
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqtt)];
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		if length(newseq(ns,:)) == 11
		seqt = newseq(ns,:);
		scorelist = [];
		allist=[];
		for t=1:length(seqt)
		seqtt = seqt;
		seqtt(t) = [];
		for s=1:length(seqtt)
			seqt3 = seqtt;
			seqt3(s) = [];
			allist = [allist; seqt3];
			scorelist= [scorelist, nwalign(seqconsensus(a_new2), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end

        a_new2=[a_new2;ali_seq];
    end
   end
  end
end

cumdist=cumsum(dist);
af=repmat('',lea,length(a_new2));
af(indexl{lmin,1},:) = a_new2(1:dist(lmin),:);
for ll=lmin+1:lmax
    af(indexl{ll,1},:)=a_new2(cumdist(ll-1)+1:cumdist(ll-1)+dist(ll),:); % final profile re-ordered
end

N = lea;
alignc = zeros(M-2,N);
for i=1:M
    for j=1:N
       alignc(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
    end
end

end
