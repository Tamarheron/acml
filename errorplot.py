def errorplot(clf, xtst):
    fig_size = [15, 6]
    plt.rcParams["figure.figsize"] = fig_size
    
    #Make a dataframe 'errors' to hold the ytest, predicted values, and errors, then sort by the error size
    errors = pd.DataFrame()
    errors['AGDIST'] = test['AGDIST']
    errors['predicted'] = clf.predict(xtst)
    errors['error'] = errors['predicted'] - errors['AGDIST']
    errors['date'] = test['strain_2_dates']
    errors = errors.sort_values('error')

    plt.figure(3)
    errors = errors.sort_values('date')
    
    #Add some random noise to the dates to avoid overplotting
    errors['jittered_dates'] = errors['date'] + np.random.uniform(-2, 2, len(test))
    
    plt.scatter(errors['jittered_dates'], errors['AGDIST'], color = 'g', s=3, marker='o')
    plt.scatter(errors['jittered_dates'], errors['predicted'], color = 'b', s=3)
    plt.scatter(errors['jittered_dates'], errors['error'], color = 'r', s=3)
    plt.xlabel('Strain 2 dates')
    plt.ylabel('Antigenic distance')
    red_patch = mpatches.Patch(color='red', label='Error')
    green_patch = mpatches.Patch(color='green', label='True value')
    blue_patch = mpatches.Patch(color='blue', label='Predicted value')
    plt.legend(handles=[red_patch, green_patch, blue_patch])
    
    plt.figure(4)
    plt.hist(errors['error'], bins=30)
    plt.xlabel('error size')
    plt.ylabel('frequency')
    plt.axvline(x=mean(errors['error']), color='r')