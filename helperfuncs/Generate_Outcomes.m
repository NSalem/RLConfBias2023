function OUT = Generate_Outcomes(ntrials,Rs,nconds)
%generates matrix of possible outcomes for good and bad option across trials
%for n conditions 
%input
%   ntrials (int) number of trials per condition, must be a multiple of 4 
%   (because probabilities used are 1/4 and 3/4)
%   Rs: (vector) two floats indicating absolute magnitude of low and high outcomes 
%output
%   OUT: matrix of size 2 x cond x trials
if nargin<2
    Rs = [0.1, 1];
end

if nargin<3
    nconds = 4;
end

%%% test preconditions
assert(numel(Rs)==2,'Rs contains wrong number of elements (should be 2)') %check we have two different possible outcomes
assert(ntrials>0 & mod(ntrials,4)==0,'number of trials is zero or not multiple of 4') %multiple of 4 is because the outcome proportions are in 1/4 - 3/4
% assert(numel(nconds)>0,'number of conditions is 0')
% nconds = 4;
lowR = Rs(1);highR = Rs(2);

%---------------------
% create trial vectors
%---------------------
OUT             = NaN(2,nconds,ntrials);
gain_A          = zeros(nconds,ntrials);
gain_B          = zeros(nconds,ntrials);
for i = 1:nconds
    if i <3 %no case >4???
        for j = 1:ntrials/4;

            gain_A(i,(j-1)*4+1:j*4) = lowR+(highR-lowR)*double(randperm(4)<4);
            gain_B(i,(j-1)*4+1:j*4) = lowR+(highR-lowR)*double(randperm(4)<2);
        end
    elseif i>2
        for j = 1:ntrials/4;
            gain_A(i,(j-1)*4+1:j*4) = -lowR-(highR-lowR)*double(randperm(4)<2);
            gain_B(i,(j-1)*4+1:j*4) = -lowR-(highR-lowR)*double(randperm(4)<4);
        end
    end
end

OUT(1,:,:)     = gain_B;%(:,randperm(size(gain_B,2))); %Better to shuffle after function
OUT(2,:,:)     = gain_A;%(:,randperm(size(gain_A,2)));

%%% test post-conditions
for icon = 1:nconds
    if icon <3
        assert(sum(OUT(2,icon,:)==highR)==ntrials*3/4 & sum(OUT(2,icon,:)==lowR)==ntrials*1/4 &...
            sum(OUT(1,icon,:)==highR)==ntrials*1/4 & sum(OUT(1,icon,:)==lowR)==ntrials*3/4,'Proportion of outcomes wrong')
    elseif icon>2
        assert(sum(OUT(2,icon,:)==-highR)==ntrials*1/4 & sum(OUT(2,icon,:)==-lowR)==ntrials*3/4&...
            sum(OUT(1,icon,:)==-highR)==ntrials*3/4 & sum(OUT(1,icon,:)==-lowR)==ntrials*1/4,'Proportion of outcomes wrong')
    end
end