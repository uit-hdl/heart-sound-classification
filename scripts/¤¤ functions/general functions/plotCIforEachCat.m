function plotCIforEachCat(cat,CI,color)
% plots confidence interval stored as column vectors in CI for each
% category in cat, where cat is an integer vector in which each integer
% represents a category, such as murmur-grade.
%% example: 
% cat = [1,2];
% CI  = [[1;2],[3.5;4]];
%% preliminary
n = length(cat);
if nargin==2
    color = strings([1,n]);
    color(:) = "k";
end
if height(CI)==3
    CI = [CI(1,:);CI(3,:)];
end
%% 
m = mean(CI);
for j=1:length(cat)
    col = color2triplet(color(j));
    line([cat(j),cat(j)],[CI(1,j),CI(2,j)],'LineWidth',2,'Color','k')
    hold on
    plot(cat(j),mean(m(j)),'o','MarkerFaceColor',col,'MarkerEdgeColor','k')
end

end