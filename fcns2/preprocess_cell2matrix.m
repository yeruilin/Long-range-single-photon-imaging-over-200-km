function [y] = preprocess_cell2matrix(MM, mm, lr, lc, bin, totDetect)

bincount=ceil((MM-mm)/10^12/bin);
y=zeros(lr,lc,bincount);

for a=1:lr
    for b=1:lc
        if size(totDetect{a,b},1)*size(totDetect{a,b},2)~=0
            y(a,b,:) = histcounts(totDetect{a,b},mm:bin*10^(12):MM);
        end
    end
end

end
