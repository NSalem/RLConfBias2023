% data informations
%% %Learning file
%% Ex: load('Sub202_Session1.mat')
% size(data) = 92 * 26  (92 trials, 26 variables)
%Column 1       2       3           4         5     
%          subj  sess   trial   con/pair  onset
%Column 6        7      8           9                     10
%         LorR  corr.   act. counter outcome confidence
%Column 11                  12     13                    14        15   
%      conf_startpoint       rt   stimA_right   stimA_Up  conf_rt
%Column 16                  17          18                    19            20   
%       Stim_left      Stim_right  act.outcome   ChosenUP  counter.outcome  
%Column 21                  22            23                    24            25  
%      UnchosenUP    choiceUP  choseA
%
%% size(stimuli) = 8 * 1  (8 symbols)
% row1: Partial Gain 75% (PC)
% row2: Partial Gain 25% (PC)
% row3: Full Gain 75% (FC)
% row4: Full Gain 25% (PC)
% row5: Partial Loss 25% (PL)
% row6: Partial Loss 75% (PL)
% row7: Full Loss 25% (FL)
% row8: Full Loss 75% (FL)
%
%%  size(money) = 1 * 4  (money from 3 sessions and total money)
%
%  Post-Learning file
% Ex: load('PostTest1.mat')
% size(data) = 92 * 26  (92 trials, 26 variables)
%Column 1       2       3                   4             5     
%          subj  trial   optionLeft   optionRight  onset
%Column 6        7             8                  9      10
%         LorR   conf   conf_startpoint      rt   conf_rt