"""
@author: mhaghigh
"""



def cross_validate(model,dataScaled,labels,groupsNum):
    """
    takes model and split object for cross validation
    
    generates various metrics for each class and each fold in a form of a dataframe
    
    """
    
    acc=[]
    predLabels=[]; #np.zeros((labels.shape[0]), dtype=bool)
    trueLabels=[];
    predProbs=[];
    trueLabelConf=[];
    i=0;
    logo = LeaveOneGroupOut()
    acc2=[]
    featureImportanceListofLists=[]
    abundance=[]
    for train_index, test_index in logo.split(dataScaled, labels, groupsNum):
    #         print("TRAIN:", train_index, "TEST:", test_index)        
        leftOutClone=groupsNum[test_index]
        X_train, X_test = dataScaled[train_index], dataScaled[test_index]        
        y_train, y_test = labels[train_index], labels[test_index]
#         sm = SMOTE()
        sm=RandomOverSampler()
#         print(i,X_train.shape,y_train.shape)

        data_res_train, labels_res_train = sm.fit_sample(X_train, y_train)
#         print(i,data_res_train.shape,labels_res_train.shape)

        model.fit(data_res_train, labels_res_train) 
        importances = model.feature_importances_
        indicesRF = np.argsort(importances)[::-1]
        topXfs=indicesRF[0:10]
#         print(np.count_nonzero(labels_res_train==True),len(labels_res_train))
        y_model = model.predict(X_test);
        prob_model=model.predict_proba(X_test)
        accuracy=accuracy_score(y_test, y_model)
        trueLabelConf0=[prob_model[i][y_test[i]] for i in range(len(y_test))];
#         accuracy=roc_auc_score(y_test, model.predict_proba(X_test))
#         accuracy=model.score(X_test,y_test)

        if accuracy>.65:
            featureImportanceListofLists.append(topXfs)
        acc.append(accuracy);
        acc2.append(str(leftOutClone[0])+':'+str(np.round(accuracy*100,2)));
#         leftOtClone
        trueLabels+=y_test.tolist();
        predLabels+=y_model.tolist();
        trueLabelConf+=trueLabelConf0;
#         predProbs+=prob_model.tolist();
#         print(trueLabelConf0,(np.sum(np.array(trueLabelConf0)>=0.95),np.sum((np.array(trueLabelConf0)>=0.95) |(np.array(trueLabelConf0)<=0.05))))
#         abnd=(np.sum(np.array(trueLabelConf0)>=0.95)/np.sum((np.array(trueLabelConf0)>=0.95) |(np.array(trueLabelConf0)<=0.05)))
#         abundance+=[abnd];
    topXfeats = set(featureImportanceListofLists[0])
    for s in featureImportanceListofLists[1:]:
        topXfeats.intersection_update(s)
    return np.mean(acc),trueLabels,predLabels,acc2,trueLabelConf,topXfeats





    
    
def balanceLeaveOneCloneOutCV(model,dataScaled,labels,groupsNum):
    """
    balance training set for each training fold, applies the classifier, evaluate roc-auc
    
    """
    from sklearn.model_selection import LeaveOneGroupOut
    from sklearn.metrics import roc_auc_score    
    
    acc=[]
    predLabels=[]; #np.zeros((labels.shape[0]), dtype=bool)
    trueLabels=[];
    predProbs=[];
    trueLabelConf=[];
    i=0;
    logo = LeaveOneGroupOut()
    acc2=[]
    featureImportanceListofLists=[]
    abundance=[]
    for train_index, test_index in logo.split(dataScaled, labels, groupsNum):
    #         print("TRAIN:", train_index, "TEST:", test_index)        
        leftOutClone=groupsNum[test_index]
        X_train, X_test = dataScaled[train_index], dataScaled[test_index]        
        y_train, y_test = labels[train_index], labels[test_index]
#         sm = SMOTE()
        sm=RandomOverSampler()
#         print(i,X_train.shape,y_train.shape)

        data_res_train, labels_res_train = sm.fit_sample(X_train, y_train)
#         print(i,data_res_train.shape,labels_res_train.shape)

        model.fit(data_res_train, labels_res_train) 
        importances = model.feature_importances_
        indicesRF = np.argsort(importances)[::-1]
        topXfs=indicesRF[0:10]
#         print(np.count_nonzero(labels_res_train==True),len(labels_res_train))
        y_model = model.predict(X_test);
        prob_model=model.predict_proba(X_test)
        accuracy=accuracy_score(y_test, y_model)
        trueLabelConf0=[prob_model[i][y_test[i]] for i in range(len(y_test))];
#         accuracy=roc_auc_score(y_test, model.predict_proba(X_test))
#         accuracy=model.score(X_test,y_test)

        if accuracy>.65:
            featureImportanceListofLists.append(topXfs)
        acc.append(accuracy);
        acc2.append(str(leftOutClone[0])+':'+str(np.round(accuracy*100,2)));
#         leftOtClone
        trueLabels+=y_test.tolist();
        predLabels+=y_model.tolist();
        trueLabelConf+=trueLabelConf0;
#         predProbs+=prob_model.tolist();
#         print(trueLabelConf0,(np.sum(np.array(trueLabelConf0)>=0.95),np.sum((np.array(trueLabelConf0)>=0.95) |(np.array(trueLabelConf0)<=0.05))))
#         abnd=(np.sum(np.array(trueLabelConf0)>=0.95)/np.sum((np.array(trueLabelConf0)>=0.95) |(np.array(trueLabelConf0)<=0.05)))
#         abundance+=[abnd];
    topXfeats = set(featureImportanceListofLists[0])
    for s in featureImportanceListofLists[1:]:
        topXfeats.intersection_update(s)
    return np.mean(acc),trueLabels,predLabels,acc2,trueLabelConf,topXfeats