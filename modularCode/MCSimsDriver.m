Initialize;
% if ~ isProdRun
%      for i=1:Ntrials
%         AttemptTrial;
%         updateTM;
%         AcceptTrial;
%         SampleTrial;
%     end
% elseif ~ usesTMMC
    for trialNo=1:Ntrials
        AttemptTrial;
        FindAccProb_L; 
        updateTM;
        AcceptTrial;
        SampleTrial;
    end
% elseif usesTMMC && ~ usesTMMCbias
%     for i=1:Ntrials
%         AttemptTrial;
%         updateTM;
%         AcceptTrial;
%         SampleTrial;
%     end
% elseif  usesTMMC && usesTMMCbias
%     for i=1:Ntrials
%         AttemptTrial;
%         updateTM;
%         AcceptTrial;
%         SampleTrial;
%     end
% end
Visualize;