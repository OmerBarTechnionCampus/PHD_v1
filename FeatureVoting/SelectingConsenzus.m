function [idx] = SelectingConsenzus(VotingMatrix,DetectionType,selectionType)
%SelectingConsenzus is a method to choose Concenzus with a variaty of
% options
%
%
%
%%
DetectionType = lower(DetectionType);
switch (DetectionType)
    case 'singlefault'
        idx = SelectingMaxima(VotingMatrix,selectionType);
    case 'multifault'
        idx = NaN;
        error('SelectingConsenzus :: selecting & detecting multiple faults Unhandeled yet')
end %swtivh type

end %SelectingConsenzus