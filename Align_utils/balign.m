function alignc = balign(seqCDR3, lf)
% note lf = N -1

Xdseq=seqCDR3;
M=length(seqCDR3);

if size(Xdseq) == [1 M],
 Xdseq = reshape(seqCDR3, [M 1]);
end

l=zeros(M,1);
lold=0;
   for i=1:M, 
     l(i)=size(Xdseq{i,1},2);   % curl brackets matrix of seqs of different length; for this object need to use cell
   end
lmax=max(l); 
lmin=min(l);

indexl=cell(lmax,1);
dist=zeros(lmax,1);
for ll=1:lmax
 indexl{ll}=find(l==ll); % indices of seqs. in each length group
 dist(ll)=length(indexl{ll}); % population of each length group
end    
             
        al=cell(lmax,1);
        p=cell(lmax,1);
        logo=cell(lmax,1);
       for ll=1:lmax 
        
        al{ll}=cell2mat(Xdseq(indexl{ll},1)); % converts a cell array into an ordinary array
        p{ll}  = seqprofile(Xdseq(indexl{ll},1),'gaps','all','counts',true); % profiles for each length
       
       end

% %create a profile for each length and then progressively align the
% %profiles old containing sequences up to lenght (L) to the new
%one containing sequences of lenght  L+1
%I start the procedure with the smallest lenght (lo)
%with a sufficient number of sequences, eg. here lo=10

nl=0; 
lo=lmin;
%lf=18; % medium-large length containing the 'bulk', i.e. the largest percentage of seqs


a_old=al{lmin,1};
%starting lenght of the alignement
l_old=dist(lmin); 
%starting profile
p_old=p{lmin};
for ll=lmin:lf
    nl=nl+1;
    if(dist(ll+1)>0)
    [pp,h1,h2]=profalign(p_old,p{ll+1}); % give optimal profile from 2 profiles
    a_new = repmat('-',l_old+dist(ll+1),ll+1);
    a_new(1:l_old,h1) = a_old;
    a_new(l_old+1:l_old+dist(ll+1),h2)=al{ll+1,1};
    else
     a_new=a_old;   
    end
    a_old= a_new;
    l_old=size(a_new,1); 
    p_old=seqprofile(a_new,'gaps','all','counts',true);    
end

hold off

%align extra sequences (longer than lf)
 %HMM Model

  Model = hmmprofstruct(lf+1);
  hmmodel=hmmprofestimate(Model, a_new); % a_new is the final alignment
 
  %Aligning other sequences to this model (removing insertions)
  a_new2=a_new;
  %ll=lf+2
  for ll=lf+2:lmax
  
    newseq=al{ll,1}; 
    for ns=1:dist(ll)
        [scores,aligned_seqs]=hmmprofalign(hmmodel,newseq(ns,:));
        %hmmodel=hmmprofestimate(Model, newseq)
        [zz,indal]=hmmprofmerge(aligned_seqs);
        ali_seq=aligned_seqs(indal);
        a_new2=[a_new2;ali_seq]; % new alignment where lengths > 18 are aligned via HMMs
    end
  end

p_new2=seqprofile(a_new2,'gaps','all','counts',true);
%reorder the alignment as the original one
%to make vector of appropriate size
cumdist=cumsum(dist);
af=repmat('',lf+1,length(a_new2));
af(indexl{lmin,1},:)=a_new2(1:dist(lmin),:);
for ll=lmin+1:lmax
    af(indexl{ll,1},:)=a_new2(cumdist(ll-1)+1:cumdist(ll-1)+dist(ll),:); % final profile re-ordered
end
    
p_f=seqprofile(af,'gaps','all','counts',true);

%Position weigth matrices with molteplicity
%alignments of zeros and ones
 N = lf+1;
        %M = size(a_new,1);
        alignc = zeros(M-2,N);
        for i=1:M
            for j=1:N
%same alphabet as matlab seqprofile used here in letter2number
               alignc(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
            end
        end

%{ 
% To write the alignment
namef = ['gen_seqs106_al_n.txt'];
fid = fopen(namef, 'a+');
for i=1:M,
 fprintf(fid, '%s \n', af(i,:));
end
fclose(fid);
%}
end
