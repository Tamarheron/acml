import numpy as np
import pandas as pd


class RegressionData:
    """
    Manage antigenic cartography data for regression analysis.

    @param filename: The C{str} name of the CSV file containing the data.
    @param dropFirstRow: Whether to drop the first row of the data
    """

    def __init__(self, filename, dropFirstRow=False, strain1ColHeading='AG1', strain2ColHeading='AG2'):
        df = pd.read_csv(filename)
        if dropFirstRow:
            self.df = df.drop(0)
        else:
            self.df = df

        #Add two columns containing 4 digit dates for each strain
        self.separate_dates(strain1ColHeading, strain2ColHeading)

        #Add two columns containing the cluster label for each strain
        self.add_clusters()

    def add_clusters(self):
        """
        Add columns that contain cluster labels, based on the 
        clusters in the map here http://www.antigenic-cartography.org/
        """
        strain_1_cluster = []
        strain_2_cluster = []

        for date in self.df['strain_1_dates']:
            if date < 1972:
                strain_1_cluster.append('HK68')
            elif date < 1975:
                    strain_1_cluster.append('EN72')
            elif date < 1977:
                    strain_1_cluster.append('VI75')
            elif date < 1979:
                    strain_1_cluster.append('TX77')
            elif date < 1987:
                    strain_1_cluster.append('BK79')
            elif date < 1989:
                    strain_1_cluster.append('SI87')
            elif date < 1992:
                    strain_1_cluster.append('BE89')
            elif date < 1995:
                    strain_1_cluster.append('BE92')
            elif date < 1997:
                    strain_1_cluster.append('WU95')
            elif date < 2002:
                    strain_1_cluster.append('SY97')
            else:
                    strain_1_cluster.append('FU02')

        for date in self.df['strain_1_dates']:
            if date < 1972:
                strain_2_cluster.append('HK68')
            elif date < 1975:
                    strain_2_cluster.append('EN72')
            elif date < 1977:
                    strain_2_cluster.append('VI75')
            elif date < 1979:
                    strain_2_cluster.append('TX77')
            elif date < 1987:
                    strain_2_cluster.append('BK79')
            elif date < 1989:
                    strain_2_cluster.append('SI87')
            elif date < 1992:
                    strain_2_cluster.append('BE89')
            elif date < 1995:
                    strain_2_cluster.append('BE92')
            elif date < 1997:
                    strain_2_cluster.append('WU95')
            elif date < 2002:
                    strain_2_cluster.append('SY97')
            else:
                    strain_2_cluster.append('FU02')

        self.df['strain_1_cluster'] = pd.Series(strain_1_cluster, index=self.df.index)
        self.df['strain_2_cluster'] = pd.Series(strain_2_cluster, index=self.df.index)

    def separate_dates(self, strain1ColHeading, strain2ColHeading):
        """
        Add columns to data frame with strain dates separated into their
        components. Adjust years from 2-digit strings to integers that
        include the century.

        @param strain1ColHeading: The C{str} strain1 label column heading in self.df.
        @param strain2ColHeading: The C{str} strain2 label column heading in self.df.
        """
        strain_1_labels = list(self.df[strain1ColHeading])
        strain_2_labels = list(self.df[strain2ColHeading])

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

    def cv(self, strain_1_cluster, strain_2_cluster):
        """
        Create a cross validation set.

        The training sets consists of all pairs where neither strain is a
        member of a particular cluster. The test sets consists of all pairs
        where at least one strain is a member of that cluster

        @params: headings of the columns which contain the labels (e.g. cluster names) that
        you want to split the data on

        @return: A cross-validation object, which is a list of train-test tuples 
        which contain the indices of the rows of the dataframe to be used in the 
        training and test data. This can be passed to scikit-learn's cross validation
        methods
        """
        df = self.df
        return [(np.where(
            (df[strain_1_cluster] != label) & (df[strain_2_cluster] != label))[0],
               np.where((df[strain_1_cluster] == label) & (df[strain_2_cluster] == label))[0])
              for label in np.unique(df[strain_1_cluster])]

    def splitByYear(self, year):
        """
        Split self.df into training and test sets by a given year.

        The training sets consist of all pairs where both strains 
        are up to the given date. The test set consists of all pairs
        where at least one strain is after the given date

        @return: A C{dict} with keys:

            'train'
            'test'
        """

        df = self.df


        train = df.loc[((df['strain_1_dates'] <= year) &
                        (df['strain_2_dates'] <= year))]

        test = df.loc[((df['strain_1_dates'] > year) |
                       (df['strain_2_dates'] > year))]

        return {

            'train': train,
            'test': test
        }
