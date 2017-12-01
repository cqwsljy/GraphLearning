function [error] = SVMforGraph(D,FD,pro)

[~,M]=size(D);


error = zeros(M,1);
for j=1:100
    index = randperm(M);
    indexTrain = index(1:round(pro*M));
    indexTest = setdiff(index,indexTrain);
    dataTrain = D(:,indexTest);
    labelTrain = FD(indexTest);
    dataTest = D(:,indexTest)
    labelTest = FD(indexTest);

    SVMStruct = fitcsvm(dataTrain',labelTrain,'Prior','uniform','Standardize',1);
    Prediciton = predict(SVMStruct,labelTest');
    error(j) = abs(labelTest- Prediciton')/length(labelTest);
end
end