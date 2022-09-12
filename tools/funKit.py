from argparse import ArgumentParser
import appConst

## Parse arguments ##
def parseArgs():
    parser = ArgumentParser()
    parser.add_argument(
        "filePath", help="File absolute path", type=str
        )
    parser.add_argument(
        "-result", "--result", help="Result file name", type=str, default='result'
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


def try1(seqs):
    self.nucInfo = []
    for i, (c1, c2) in enumerate(zip(*seqs), 1):
        self.nucInfo.append(
            {
                'index' : i,
                'ref' : c1,
                'que' : c2
            }
        )
    print()
    print("Count of Mismatches is :", len(mismatch))
    print("Count of N is :", len(ns))
    print("Position of Mismatches :", *mismatch)
    print("Position of N :", *ns)
