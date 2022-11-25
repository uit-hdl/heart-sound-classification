function class = num2class(x)
if  x==-1
    class = "normal";
elseif x==0
    class = "uncertain";
else
    class = "abnormal";
end
end