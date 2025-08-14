function out = TwinFrequency2(input,maxT,err)
%input a list of misorientations and the error allowed in deltaA to count
%as a twin, output frequencies of each twin
n = size(input,1);
twins = zeros([n 1]);
twins(input < 63.8+err & input > 63.8-err) = 1;
if maxT >= 2
    twins(input < 52.4+err & input > 52.4-err) = 2;
end
if maxT >= 3
    twins(input < 11.4+err & input > 11.4-err) = 3;
end
if maxT >= 4
    twins(input < 75.2+err & input > 75.2-err) = 4;
end
if maxT >= 5
    twins(input < 41+err & input > 41-err) = 5;
end
if maxT >= 6
    twins(input < 22.8+err & input > 22.8-err) = 6;
end
if maxT >= 7
    twins(input < 86.6+err & input > 86.6-err) = 7;
end
if maxT >= 8
    twins(input < 29.6+err & input > 29.6-err) = 8;
end
if maxT >= 9
    twins(input < 34.2+err & input > 34.2-err) = 9;
end
if maxT >= 10
    twins(input < 82+err & input > 82-err) = 10;
end

out = zeros([maxT+1 1]);
for a = 1:maxT+1
    out(a) = sum(twins == a-1);
end