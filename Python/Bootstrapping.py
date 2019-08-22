from itertools import repeat
from numpy import array, append, empty, std, mean
from numpy.random import rand, choice
from collections import Counter
from scipy.stats import sem, t
from pandas import DataFrame, read_excel

def ConfidenceInterval(data, confidence=0.95):
    a = 1.0 * array(data)
    n = len(a)
    se = sem(a)
    h = se * t.ppf((1 + confidence) / 2., n-1)
    return h

def bootstrap_resample(X, n=None):
    """ Bootstrap resample of array_like data
    
    Parameters
    ----------
    X : array_like
      data to resample
    n : int, optional
      length of resampled array, equal to len(X) if n==None
      
    Results
    -------
    returns X_resamples
    """
    
    if n == None:
        n = len(X)
        
    resample_i = choice(X,size=n, replace=True)
    return resample_i

Sample = []
MutatedPositions = ["11", "12", "13", "15", "16", "17", "18", "34", "35", "36", "37", "39"]
OutsideMutatedPositionsName = "Outside Mutated Positions"
DoubleMutantsName = "Double Mutants"
filename = r"D:\Promotion\Matlab\Data\NGS_6.2Mix_sorted_Ch_aMT\2Mix\Combined\2Mix_singles_Si.xlsx"

#Import Excel File and convert into a dictionary
DataExcel = read_excel(filename) 
SampleDataFrame = DataFrame(DataExcel, columns= ["Mutant","Number of mutant in sorted population"])


UnsortedUniqueElements = array(SampleDataFrame["Mutant"])
ElementCount = array(SampleDataFrame["Number of mutant in sorted population"])
NumberOfDoubleMutants = read_excel(filename, 'Sheet1', usecols = "F", nrows=1, header =None).iat[0,0] - sum(ElementCount)

OutsideMutatedPositions = 0
for i in range(len(ElementCount)):
    if UnsortedUniqueElements[i][1:-1] in MutatedPositions or UnsortedUniqueElements[i] == "WT":
        Sample.extend(repeat(UnsortedUniqueElements[i], ElementCount[i]))
    else:
        OutsideMutatedPositions += ElementCount[i]
Sample.extend(repeat(OutsideMutatedPositionsName,OutsideMutatedPositions))
Sample.extend(repeat(DoubleMutantsName,NumberOfDoubleMutants))
    
#Sorting Mutants
SampleCounter = Counter(Sample)
SortedUniqueElements = list(SampleCounter.keys())
SortedUniqueElements.remove("WT")
SortedUniqueElements.remove(OutsideMutatedPositionsName)
SortedUniqueElements.remove(DoubleMutantsName)
SortedUniqueElements.sort(key=lambda x: (int(x[1:-1]), x[-1]))
SortedUniqueElements.append(OutsideMutatedPositionsName)
SortedUniqueElements.append(DoubleMutantsName)
SortedUniqueElements.insert(0,"WT")

#Bootstrapping
AllResamplesCounter = empty([0,len(SortedUniqueElements)])
for _ in range(1000):
    Resample = bootstrap_resample(array(Sample))
    ResampleCounter = Counter(Resample)
    SortedResample = array([])
    for i in SortedUniqueElements:
        #if i not in ResampleCounter:
           # ResampleCounter[i] = 0
        SortedResample = append(SortedResample,ResampleCounter[i]/len(Resample))
    AllResamplesCounter = append(AllResamplesCounter,[SortedResample],axis=0)

#Calculating Uncertainties
DictAllIntervals = {}
DictStandardDeviations = {}
DictMeans = {}
DictStandardDeviationsAndConfidenceIntervals = {}
for i in range(len(SortedUniqueElements)):
    DictAllIntervals[SortedUniqueElements[i]] = [std(AllResamplesCounter[:,i]),sem(AllResamplesCounter[:,i]),ConfidenceInterval(AllResamplesCounter[:,i]),mean(AllResamplesCounter[:,i]),SampleCounter[SortedUniqueElements[i]]/len(Sample)]
    DictStandardDeviations[SortedUniqueElements[i]] = std(AllResamplesCounter[:,i])
    DictStandardDeviationsAndConfidenceIntervals[SortedUniqueElements[i]] = [std(AllResamplesCounter[:,i]),ConfidenceInterval(AllResamplesCounter[:,i])]
    DictMeans[SortedUniqueElements[i]] = mean(AllResamplesCounter[:,i])

#Export Results to Excel file
DictToExport = DictAllIntervals
ResultsDataFrame = DataFrame(DictToExport.values(), columns= ["Standard Deviation","Standard Error","95% Confidence Interval","Mean Frequency","Measured Frequency"], index = DictToExport.keys())
ResultsDataFrame.to_excel(r"D:\Promotion\Bootstrapping\2Mix_presorted_BovineTrypsin.xlsx")

print ("Calculation is finished")