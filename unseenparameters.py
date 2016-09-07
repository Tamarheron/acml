def unseenParameters(clf, xtrain):
    unseen = []
    for column in xtrain.columns:
        if sum(xtrain[column])==0:
            unseen.append(column)
    params = pd.DataFrame(data = clf.coef_, index = xtrain.columns, columns=('Coefficient',))
    return(params.loc[unseen])