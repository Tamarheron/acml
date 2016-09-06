import numpy as np
import pandas as pd


class RegressionData:
    """
    Manage antigenic cartography data for regression analysis.

    @param filename: The C{str} name of the CSV file containing the data.
    """

    def __init__(self, filename, dropFirstRow=False):
        df = pd.read_csv(filename)
        if dropFirstRow:
            self.df = df.drop(0)
        else:
            self.df = df

        self.separate_dates('AG1', 'AG2')
        self.add_clusters('strain_1_dates', 'strain_2_dates')

    def add_clusters(self, date1, date2):
        """
        Add clusters....
        """
        cluster1 = []
        cluster2 = []

        for date in self.df[date1]:
            if date < 1972:
                cluster1.append('HK68')
            elif date < 1975:
                    cluster1.append('EN72')
            elif date < 1977:
                    cluster1.append('VI75')
            elif date < 1979:
                    cluster1.append('TX77')
            elif date < 1987:
                    cluster1.append('BK79')
            elif date < 1989:
                    cluster1.append('SI87')
            elif date < 1992:
                    cluster1.append('BE89')
            elif date < 1995:
                    cluster1.append('BE92')
            elif date < 1997:
                    cluster1.append('WU95')
            elif date < 2002:
                    cluster1.append('SY97')
            else:
                    cluster1.append('FU02')

        for date in self.df[date2]:
            if date < 1972:
                cluster2.append('HK68')
            elif date < 1975:
                    cluster2.append('EN72')
            elif date < 1977:
                    cluster2.append('VI75')
            elif date < 1979:
                    cluster2.append('TX77')
            elif date < 1987:
                    cluster2.append('BK79')
            elif date < 1989:
                    cluster2.append('SI87')
            elif date < 1992:
                    cluster2.append('BE89')
            elif date < 1995:
                    cluster2.append('BE92')
            elif date < 1997:
                    cluster2.append('WU95')
            elif date < 2002:
                    cluster2.append('SY97')
            else:
                    cluster2.append('FU02')

        self.df['cluster1'] = pd.Series(cluster1, index=self.df.index)
        self.df['cluster2'] = pd.Series(cluster2, index=self.df.index)

    def separate_dates(self, col1, col2):
        """
        Add columns to data frame with strain dates separated into their
        components. Adjust years from 2-digit strings to integers that
        include the century.

        @param col1: The C{str} strain1 label column heading in self.df.
        @param col2: The C{str} strain2 label column heading in self.df.
        """
        strain_1_labels = list(self.df[col1])
        strain_2_labels = list(self.df[col2])

        strain_1_dates = []
        strain_2_dates = []

        for label in strain_1_labels:
            date = int(list(label.split('/'))[-1])
            if date < 60:
                date += 2000
            else:
                date += 1900
            strain_1_dates.append(date)

        for label2 in strain_2_labels:
            date = int(list(label2.split('/'))[-1])
            if date < 60:
                date += 2000
            else:
                date += 1900
            strain_2_dates.append(date)

        self.df['strain_1_dates'] = pd.Series(strain_1_dates,
                                              index=self.df.index)
        self.df['strain_2_dates'] = pd.Series(strain_2_dates,
                                              index=self.df.index)

    def cv(self, cluster1, cluster2):
        """
        Create a cross validation set.

        The training set consists of all pairs where neither strain is a
        member of a particular cluster. The test set consists of all pairs
        where at least one strain is a member of that cluster

        @return: A cross-validation something...
        """
        df = self.df
        return [(np.where(
            (df[cluster1] != label) & (df[cluster2] != label))[0],
               np.where((df[cluster1] == label) & (df[cluster2] == label))[0])
              for label in np.unique(df[cluster1])]

    def splitByYear(self, year):
        """
        Split self.df by a given year.

        @return: A C{dict} with keys:

            'X': the X data...
        """

        df = self.df
        X = df.drop(['AG1', 'AG2', 'AG-DIST', 'NUM-MUTATIONS',
                     'strain_1_dates', 'strain_2_dates', 'cluster1',
                     'cluster2'], axis=1)

        y = df['AG-DIST']

        train = df.loc[((df['strain_1_dates'] <= year) &
                        (df['strain_2_dates'] <= year))]

        test = df.loc[((df['strain_1_dates'] > year) |
                       (df['strain_2_dates'] > year))]

        xtrain = train.drop(
            ['AG1', 'AG2', 'AG-DIST', 'NUM-MUTATIONS',
             'strain_1_dates', 'strain_2_dates', 'cluster1', 'cluster2'],
            axis=1)

        ytrain = train['AG-DIST']

        xtest = test.drop(
            ['AG1', 'AG2', 'AG-DIST', 'NUM-MUTATIONS',
             'strain_1_dates', 'strain_2_dates', 'cluster1', 'cluster2'],
            axis=1)

        ytest = test['AG-DIST']

        return {
            'X': X,
            'y': y,
            'train': train,
            'test': test,
            'xtrain': xtrain,
            'ytrain': ytrain,
            'xtest': xtest,
            'ytest': ytest,
        }
