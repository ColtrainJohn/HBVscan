from argparse import ArgumentParser
import pandas as pd



## HBVisor
from . import appConst

## Parse arguments ##
def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "data", help="File absolute path", type=str
    )
    parser.add_argument(
        "out_all", help="File absolute path", type=str, default=False
    )
    parser.add_argument(
        "out_report", help="File absolute path", type=str, default=False
    )
    args = parser.parse_args()
    return args


def parseBlastOut(blastOut):
    blastFeatures = {}
    for key, value in zip(
            *[
        appConst.blastOutCols.split(' '),
        blastOut.split('\t')
            ]):
        blastFeatures.update({key : value})
    return blastFeatures


def checkPositions(resultDict):
    nucInfo = []
    for i, (c1, c2) in enumerate(zip(
        *[
            resultDict['sseq'],
            resultDict['qseq']
                                ]), 0):
        nucInfo.append(
            {
                'indexRef' :  int(resultDict['sstart']) - 1 + i,
                'indexQuery' : int(resultDict['qstart']) - 1 + i,
                'ref' : c1,
                'que' : c2
            }
        )
    tab = pd.DataFrame(nucInfo)
    tab = tab.loc[tab['ref'] != tab['que']]
    tab = tab.reset_index()[['indexRef', 'indexQuery', 'ref', 'que']]
    return tab

