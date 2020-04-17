function [pwmat] = seqprofile_rew(af,q,We)


[M,N]=size(af);
msa21=zeros(M,N*q);

align=zeros(M,N);
for i=1:M
    for j=1:N
       align(i,j) = letter2number(af(i,j)); % convert into numbers to count them and obtain frequencies
end
end

    for i=1:N;
        for a=1:q; % a.a. are labeled from 1 to q
	msa21(:,(i-1)*q+a)=(align(:,i)==a);
        end
    end

Meff = sum(We);
w = We'/Meff;
if size(msa21,1) == 1
 pw = w.*msa21;
else
pw = sum(w.*msa21);
end
pwmat = reshape(pw,q,N);
end

