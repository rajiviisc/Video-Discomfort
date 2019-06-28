function selective_scores = calculate_mos(fName,pathVar,sublist)

load(strcat(fName,pathVar,'zscores.mat'));
scores = zeros(size(zscoreMat,2),1);
count = zeros(size(zscoreMat,2),1);
for i=1:size(sublist,2)  

    idx = find(~isnan(zscoreMat(sublist(i),:)));
    z_score = zscoreMat(sublist(i),idx);

    for j=1:size(z_score,2)
        scores(idx(j)) = scores(idx(j)) + z_score(j);
        count(idx(j)) = count(idx(j))+1;        
    end    
end
selective_scores = scores./count; 

