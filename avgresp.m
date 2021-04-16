function resp_coord = avgresp

numStim = 4; % 5 including image stimuli
contlist = [2 8 16];
resp_coord = zeros(numStim,3,10);
ratings = zeros(numStim,3,10,40);
ratingIndex = ones(numStim,3,10);
for i = 1:10       % Loop through all 10 subjects
    name = sprintf('respMat_s%d.mat',i);
    currentStruct = load(name);         % Load response, but the cell is loaded as a struct
    currentMat = currentStruct.respMat; % Convert the loaded struct back to cell
    for k = 1:520
        currentStim = currentMat{2}.stimuli(k);
        currentCont = currentMat{2}.contrast(k);
        currentRating = currentMat{2}.rate_percept(k);
        
        for stim = 1:numStim
            
            if stim == currentStim
            for cont = 1:3
                currentContIndex = find(currentCont == contlist);
                if cont == currentContIndex
                    ratings(stim,cont,i,ratingIndex) = currentRating;
                    ratingIndex(stim,cont,i) = ratingIndex(stim,cont,i)+1;
                end
            end
            end
        end
    end
end

for Stim = 1:numStim
    for Cont = 1:3
        for Sub = 1:10
            resp_coord(Stim,Cont,Sub) = mean(ratings(Stim,Cont,Sub,:));
        end
    end
end
end

    
                    
        