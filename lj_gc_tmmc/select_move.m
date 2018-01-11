function [mSelect] = select_move(moveProb)
mTypes=length(moveProb);
selRand = rand();
for i = 1:mTypes
    if selRand<moveProb(i)
        mSelect = i;
        return;
    end
end

end