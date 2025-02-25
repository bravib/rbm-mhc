function alignc_new = balign_peps_seed(seqs_seed, lea, lepmin, lepmax, seqs_new, weights)


%%%% Build HMM profile from seqs_seed %%%%
yes_weight = 0;
if nargin > 5 & weights
 weights = importdata('weights.txt');
 Weights = weights';
 yes_weight = 1;
end

yes_red = 0;
sca=log(2);
%sca=1;
SM = 'BLOSUM62'; % default BLOSUM50
GO = 8; % default 8

Xdseq = seqs_seed;
M=length(seqs_seed);
if size(Xdseq) == [1 M],
 Xdseq = reshape(seqs_seed, [M 1]);
end

l=zeros(M,1);
lold=0;
   for i=1:M, 
     l(i)=size(Xdseq{i,1}, 2); 
 end
lmax=max(l); 
lmin=min(l);

if lmin==lmax,
	af = cell2mat(Xdseq);
	N = length(af(1,:));
	alignc = zeros(M,N);
	for i=1:M
	    for j=1:N
	       alignc(i,j) = letter2number(af(i,j));
	    end
	end

alignc_seed = alignc;
af_seed = af;
a_new = af;
a_new2= af;
Model = hmmprofstruct(lmin);
hmmodel = hmmprofestimate(Model, a_new);
a_new_repr = a_new2;
else
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
p{ll} = seqprofile_rew(al{ll}, 21, W{ll});
else
p{ll} = seqprofile(Xdseq(indexl{ll},1),'gaps','all','counts',true); % seq profiles (i.e. aa count in each position) for each length
end
end

if lepmax > lmax
lepmax = lmax;
end
if lepmin < lmin
lepmin = lmin;
end

nl=0; 
lo=lepmin; 
a_old=al{lepmin,1};
a_new2=a_old;
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
	    [pp,h1,h2] = profalign(p_old, p{ll}, 'ScoringMatrix', SM);
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
else
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
a_new2= a_old;
end

Model = hmmprofstruct(lea);
hmmodel = hmmprofestimate(Model, a_new);
a_new_repr = a_new2; % representative 9-long subalignment used in HMM and here for alignment scores (with 8 and 9)

if lea < lmax
if lea + 1 < lmax
  for ll = lea + 1 : lmax
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
 else
    ll=lmax;
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
if lmin < lepmin
if lmin < lepmin - 1
  for ll = lmin:lepmin-1
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
 else
    ll=lmin;
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

cumdist=cumsum(dist);
af=repmat('',lea,length(a_new2));
af(indexl{lmin,1},:) = a_new2(1:dist(lmin),:);
for ll=lmin+1:lmax
    af(indexl{ll,1},:) = a_new2(cumdist(ll-1)+1:cumdist(ll-1)+dist(ll),:); % final profile re-ordered
end

N = lea;
alignc = zeros(M-2,N);
for i=1:M
    for j=1:N
       alignc(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
    end
end
alignc_seed = alignc;
af_seed = af;

end

%a_new_repr = al{lea};
%%%%%% From here align new sequences to the same seed profile %%%%%%

Xdseq = seqs_new;
M=length(seqs_new);
if size(Xdseq) == [1 M],
 Xdseq = reshape(seqs_new, [M 1]);
end
l=zeros(M,1);
lold=0;
   for i=1:M, 
     l(i)=size(Xdseq{i,1}, 2);   % Length distribution of new sequences
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
alignc_new = alignc;
af_new = af;

scores_new = ones(M,1);
namef=[rootf,'/Align_utils/scores.txt'];
fid = fopen(namef, 'w');
for i=1:length(scores_new)
 fprintf(fid,'%f\n', scores_new(i));
end
fclose(fid);

scores_new_cons= ones(M,1);
namef=[rootf,'/Align_utils/scores_cons.txt'];
fid = fopen(namef, 'w');
for i=1:length(scores_new_cons)
 fprintf(fid,'%f\n', scores_new_cons(i));
end
fclose(fid);

return
end


indexl=cell(lmax,1);
dist=zeros(lmax,1);
for ll=1:lmax
 indexl{ll}=find(l==ll); % indices of seqs. in each length group
 dist(ll)=length(indexl{ll}); % population of each length group
end

al=cell(lmax,1);
p=cell(lmax,1);
for ll=1:lmax 
al{ll} = cell2mat(Xdseq(indexl{ll},1)); % converts a cell array into an ordinary array
p{ll}  = seqprofile(Xdseq(indexl{ll},1),'gaps','all','counts',true); % seq profiles (i.e. aa count in each position) for each length
end


% align new sequences to the HMM model
LA = length(a_new2);
a_new3 = a_new2;
%a_new_repr = a_new2; % probably the most consistent choice is to align to the full final sample?

scores_new0 = [];
scores_new0_cons = [];
if lmin < lmax % 'align', using the HMM probabilistic model, only if the sequences differ from the average length of the profile

for ll = lmin:lmax
newseq = al{ll,1}; 
 if ll ~= lea
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
		 scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
			scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end

		a_new3=[a_new3;ali_seq]; % new alignment
		scores_new0 = [scores_new0, scores];
		%scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), newseq(ns,:), 'Scale', log(2))*length(newseq(ns,:))];
		scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), newseq(ns,:), 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];
%{		
		rf=[];
		for u = 1:size(a_new_repr,1),
			rf = [rf, nwalign(a_new_repr(u,:), newseq(ns,:), 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];
		end
		scores_new0_cons = [scores_new0_cons, max(rf)];
%}
	end
else

for ns=1:dist(ll)
	   ali_seq=newseq(ns,:);

	   if length(newseq(ns,:)) ~= 9 & length(strfind(ali_seq,'-')) > 1
		if length(newseq(ns,:)) == 8
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)-1
		 seqtt = [seqt(1:t), '-', seqt(t+1:end)];
		 allist= [allist; seqtt];
		 scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
			scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end
  		a_new3=[a_new3;ali_seq]; % new alignment
		[scores,~]=hmmprofalign(hmmodel,ali_seq);
		scores_new0 = [scores_new0, scores];
		scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), ali_seq, 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];

		end
 end
end

else
ll = lmin;

if ll ~= lea
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
		 scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
			scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end

		a_new3=[a_new3;ali_seq]; % new alignment
		scores_new0 = [scores_new0, scores];
		%disp(seqconsensus(a_new2))
		scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), newseq(ns,:),'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];
	end
else
for ns=1:dist(ll)
	   ali_seq=newseq(ns,:);
	   if length(newseq(ns,:)) ~= 9 & length(strfind(ali_seq,'-')) > 1 & yes_red
		if length(newseq(ns,:)) == 8
		seqt = newseq(ns,:);
		scorelist = [];
		allist = [];
		for t=1:length(seqt)-1
		 seqtt = [seqt(1:t), '-', seqt(t+1:end)];
		 allist= [allist; seqtt];
		 scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		end
		seqtt = ['-', seqt];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
		seqtt = [seqt, '-'];
		allist= [allist; seqtt];
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
		scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqtt)];
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
			scorelist= [scorelist, nwalign(seqconsensus(a_new3), seqt3)];
		 end
		end
		[m,p] = max(scorelist);
		ali_seq = allist(p,:);
		end

		end
	   a_new3=[a_new3;ali_seq]; % new alignment
	   [scores,~]=hmmprofalign(hmmodel,ali_seq);
	     scores_new0 = [scores_new0, scores];
	     scores_new0_cons = [scores_new0_cons, nwalign(seqconsensus(a_new_repr), ali_seq, 'Scale', sca, 'ScoringMatrix', SM, 'GapOpen', GO)];

	end
end
end
a_new4 = a_new3((LA + 1) : (LA + length(seqs_new)),:);
cumdist = cumsum(dist);
af = repmat('', lea, size(a_new4,1));
scores_new = zeros(size(a_new4,1),1);
scores_new_cons = zeros(size(a_new4,1),1);
af(indexl{lmin,1},:)=a_new4(1:dist(lmin),:);
scores_new(indexl{lmin,1}) = scores_new0(1:dist(lmin));
scores_new_cons(indexl{lmin,1}) = scores_new0_cons(1:dist(lmin));
for ll=lmin+1:lmax
   af(indexl{ll,1},:) = a_new4(cumdist(ll-1)+1:cumdist(ll-1)+dist(ll),:); % final profile re-ordered
   scores_new(indexl{ll,1}) = scores_new0(cumdist(ll-1)+1:cumdist(ll-1) + dist(ll)); % final HMM optimal scores re-ordered
   scores_new_cons(indexl{ll,1}) = scores_new0_cons(cumdist(ll-1)+1:cumdist(ll-1) + dist(ll)); % final alignment scores
end
N = lea;
alignc = zeros(M-2,N);
for i=1:M
    for j=1:N
       alignc(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
    end
end
alignc_new = alignc;
af_new = af;

namef=[rootf,'/Align_utils/scores.txt'];
fid = fopen(namef, 'w');
for i=1:length(scores_new)
 fprintf(fid,'%f\n', scores_new(i));
end
fclose(fid);
namef=[rootf,'/Align_utils/scores_cons.txt'];
fid = fopen(namef, 'w');
for i=1:length(scores_new_cons)
 fprintf(fid,'%f\n', scores_new_cons(i));
end
fclose(fid);

end
