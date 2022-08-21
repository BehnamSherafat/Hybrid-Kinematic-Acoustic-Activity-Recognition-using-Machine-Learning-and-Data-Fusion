load fisheriris

rng('default') % For reproducibility
c = cvpartition(species,'KFold',5);

discrCVModel = fitcdiscr(meas,species,'CVPartition',c);
treeCVModel = fitctree(meas,species,'CVPartition',c);

discrRate = kfoldLoss(discrCVModel)
discrRate = 0.0200
treeRate = kfoldLoss(treeCVModel)