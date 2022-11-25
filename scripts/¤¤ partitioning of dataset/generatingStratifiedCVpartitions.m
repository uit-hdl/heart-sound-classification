% In this script, I generate the cross-validation splits which is used in
% all other code. For this, I use the matlab official function cvpartition.
% I define categories that I wish to balance across, such as disease cases
% and presence of murmur. Then, generate a large number of CV-partitions,
% and select the one with the most even distribution (in an average sense)
% across the classes. The result is saved in the structure "CVpartitions"
% found in the "saved variables" folder.
%%
% note that HSdataTrain = HSdata(union(Jtrain0,Jval0),:);
data0 = HSdataTrain;
data0 = data0(randperm(height(data0)),:);

Jall = 1:height(data0);
% classes to balance across CV-partitions:
Jmur1 = data0.Murmur_1_grade_ref_ny_T72>=1;
Jmur2 = data0.Murmur_2_grade_ref_ny_T72>=1;
Jmur3 = data0.Murmur_3_grade_ref_ny_T72>=1;
Jmur4 = data0.Murmur_4_grade_ref_ny_T72>=1;
Jas   = data0.ASPRESENCE_T72==1;
Jms   = data0.MSPRESENCE_T72==1;
Jar   = data0.ARGRADE_T72>=3;
Jmr   = data0.MRGRADE_T72>=3;

% folds = cell(8,8);
% c = cvpartition(Jms,'KFold',8);
% for i=1:8
%     folds{1,i} = find(data0.MSPRESENCE_T72(c.test(i))==1);
% end


minUnequalness = 2;
warning('off')
% generate a large number of cv-partitions, and keep the one which has the
% most balanced class-proportions.
for j=1:5000000
% c = cvpartition(sum([Jmur1';Jmur2';Jmur3';Jmur4';Jas';Jms';Jar';Jmr']),'KFold',8);
% generate random cv-partition (balanced across AS cases)
c = cvpartition(Jas,'KFold',8);
% sum(data0.ASPRESENCE_T72(c.test(4)))

p = zeros(5,8);
% compute class-representations in each partition:
for i=1:8
    p(1,i) = sum(data0.ARGRADE_T72(c.test(i))>=3);
    p(2,i) = sum(data0.MRGRADE_T72(c.test(i))>=3);
    p(3,i) = sum(data0.ASGRADE_T72(c.test(i))>=1);
    p(4,i) = sum(data0.MSGRADE_T72(c.test(i))>=1);
    p(5,i) = sum(data0.Murmur_1_grade_ref_ny_T72(c.test(i))>=1) + ...
                 sum(data0.Murmur_2_grade_ref_ny_T72(c.test(i))>=1) + ...
                 sum(data0.Murmur_3_grade_ref_ny_T72(c.test(i))>=1) + ...
                 sum(data0.Murmur_4_grade_ref_ny_T72(c.test(i))>=1);
end

spread = range(p')./mean(p');

% measure the eveness of the spread:
unequalness = .1*spread(1)+.1*spread(2)+.4*spread(3)+.25*spread(4)+.15*spread(5);

if unequalness<minUnequalness
    minUnequalness = unequalness;
    bestP = p  %#ok<*NOPTS>
    bestPartition = c;
end

end


bestP
    
clf
for i=1:5
    subplot(3,2,i);
    bar(bestP(i,:))
end
range(bestP')./min(bestP') %#ok<*UDIM>

%% save the most balanced CV-partition 
% CVsplotBal.train = cell(8,1);
% CVsplotBal.val   = cell(8,1);
for i=1:8
    CVsplitBal.train.id{i,1} = ind2id(bestPartition.training(i),data0);
    CVsplitBal.train.I{i,1}  = id2ind(CVsplitBal.train.id{i,1},HSdata);
    CVsplitBal.train.J{i,1}  = find(CVsplitBal.train.I{i,1});
    
    CVsplitBal.val.id{i,1}  = ind2id(bestPartition.test(i),data0);
    CVsplitBal.val.I{i,1}   = id2ind(CVsplitBal.val.id{i,1},HSdata);
    CVsplitBal.val.J{i,1}   = find(CVsplitBal.val.I{i,1});
end

