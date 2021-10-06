function [scores] = sort_by_choice_order(matrix, choiceOrder)

    nSubject = size(matrix,1);
    nChoice = size(matrix,2);
    scores = zeros(nSubject, nChoice);
    
    for i = 1:nSubject
        sc = matrix(i,:);
        ch = choiceOrder(i,:);
        [~ , sortedChoiceOrder] = sort(ch);  %sort choice order from 1 to 6
        scores(i,:) = sc(ch(sortedChoiceOrder));    %rearrange scores by sorted choice order
    end
    
    scores(isnan(scores))=0;
    
end

