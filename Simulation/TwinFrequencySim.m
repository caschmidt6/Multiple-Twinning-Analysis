function out = TwinFrequencySim(input)

out = zeros(5,1);

[hc,~] = histcounts(input,0:1:90);

out(2,1) = sum(hc(1,65));
out(3,1) = sum(hc(1,53));
out(4,1) = sum(hc(1,13));
out(5,1) = sum(hc(1,77));
out(1,1) = sum(hc) - sum(out);