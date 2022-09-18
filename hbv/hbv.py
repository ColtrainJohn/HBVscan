# default
from subprocess import run
from argparse import ArgumentParser
import sys, pathlib #!
import itertools


# special
from Bio import SeqIO
import pandas as pd


# HBVisor
from . import appConst
from . import funKit



## Custom SEQ class ##
class newHBVpart:
    def __init__(self, inFile, outAll, outReport):
        self.filePath = inFile
        self.root = str(pathlib.Path(__file__).parent.resolve())
        self.bioSeqIO = SeqIO.read(inFile, 'fasta')
        self.blastOn()
        self.posTab = funKit.checkPositions(self.blastFeatures)
        self.reportTab = self.makeReport()
        ## Output files
        if not (outAll and outReport):
            print(self.posTab)
            print(self.reportTab)
            print(dict((k, self.reportDict[k]) for k in ['Name', 'Conclusion'] if k in self.reportDict))
        if outAll:
            self.posTab.to_csv(outAll + '.csv')
        if outReport:
            self.reportTab.to_csv(outReport + '.csv')


    def blastOn(self):
        appConst.argsCmd.update({'-db' : appConst.argsCmd['-db'].substitute(root=f'{self.root}' + '/hbv1')})
        appConst.argsCmd.update({'-query' : appConst.argsCmd['-query'].substitute(file=self.filePath)})
        self.blastDone = run(
            [
                appConst.blastCmd.substitute(root=f'{self.root}' + '/blastn'),
                    *list(itertools.chain.from_iterable(
                        zip(
                            list(appConst.argsCmd.keys()), 
                            list(appConst.argsCmd.values())
                        )
                    )
                )
            , '-subject_besthit'],
            capture_output=True
        ).stdout.decode('utf-8').strip('\n')
        self.blastFeatures = funKit.parseBlastOut(self.blastDone)

    def makeReport(self):
        self.reportDict = {
                'Name' : self.bioSeqIO.description,
                'Conclusion' : self.blastFeatures['sseqid'],
                'refStart' : self.blastFeatures['sstart'],
                'refEnd' : self.blastFeatures['send'],
                'queStart' : self.blastFeatures['qstart'],
                'queEnd' : self.blastFeatures['qend'],
                'mismatchCounts' : self.blastFeatures['mismatch'],
                'gapCounts' : self.blastFeatures['gaps'].strip('\n')
            }
        tab = pd.DataFrame.from_dict(self.reportDict, 'index')
        tab.columns = ['Result']
        return tab

    def checkPositions(self):
        self.nucInfo = []
        for i, (c1, c2) in enumerate(zip(
            *[
                self.blastFeatures['sseq'], 
                self.blastFeatures['qseq']
                                            ]), 0):
            self.nucInfo.append(
                {
                    'indexRef' :  int(self.blastFeatures['sstart']) - 1 + i,
                    'indexQuery' : int(self.blastFeatures['qstart']) - 1 + i,
                    'ref' : c1,
                    'que' : c2
                                    } 
            )
        tab = pd.DataFrame(self.nucInfo)
        tab = tab.loc[tab['ref'] != tab['que']]
        tab = tab.reset_index()[['indexRef', 'indexQuery', 'ref', 'que']]
        return tab


if __name__ == "__main__":
    args = funKit.parseArgs()
    hand = newHBVpart(args.inFile, args.outAll, args.outReport)
