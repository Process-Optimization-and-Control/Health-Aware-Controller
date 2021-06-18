%test

x1 = ones(9,1);

[y1c] = NNFit(x1)

an= [];
for ii = 1:9
    temp = zeros(9,1);
    
    temp(ii) = 1;
    
    [y1] = NNFit(x1 + temp*0.05);
    
    an = [an; (y1 - y1c)/0.05];
    
end

an