from unittest import TestCase

# import pandas as pd

from regressiondata import RegressionData


class TestRegressionDataForUpTo10Mutations(TestCase):

    def setUp(self):
        """
        Read the CSV file.
        """
        DATAFILE = ('antigenic-genetic-pret8-positional-no-0-and-upto-'
                    '10-mutations.csv')
        self.data = RegressionData(DATAFILE, dropFirstRow=True)

    def testShape(self):
        """
        The df variable must have the expected shape.
        """
        self.assertEqual((6723, 307), self.data.df.shape)

    def testLength(self):
        """
        The df variable must have the expected length.
        """
        self.assertEqual(6723, len(self.data.df))

    def testSplitByYearKeys(self):
        """
        The splitByYear function must return a dictionary with the expected
        keys.
        """
        expected = sorted(['test', 'xtest', 'ytrain', 'ytest', 'y', 'X',
                           'train', 'xtrain'])

        self.assertEqual(expected, sorted(self.data.splitByYear(1994).keys()))

    def testSplitByYear_X(self):
        """
        The splitByYear function must return an X item with 299 columns when
        splitting on 1994.
        """
        self.maxDiff = None
        self.assertEqual((6723, 299),
                         self.data.splitByYear(1994)['X'].shape)


# class TestSeparateDates(TestCase):

#     def testNonExistentColumn(self):
#         """
#         The separate_dates function must raise a KeyError if it is passed
#         a non-existent column name.
#         """
#         dataFrame = pd.DataFrame({}, columns=[])
#         self.assertRaises(KeyError, separate_dates, dataFrame, 'x', 'y')
