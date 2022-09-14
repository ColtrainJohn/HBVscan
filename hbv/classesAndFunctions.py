import operator
from Bio import SeqIO
from Bio import pairwise2

#class Classificator:
#    def __init__(self, refsPath):


class hbvRefs:
    def __init__(self, refsPath):
        self.refs = list(SeqIO.parse(refsPath, 'fasta'))
        self.genoDict = dict()
        
        for ref in self.refs:
            self.genoDict.update(
                    {
                        ref.description : str(ref.seq)
                    }
                                )
    def throwOn(self, fasta):
        query = SeqIO.parse(fasta, 'fasta')
        outDict = dict()
        for each in query:
            print('Query name: ' + each.id)
            maxScore = 1
            finalAln = '?'
            gen = '?'
            for refItem in self.genoDict.items():
                aln = pairwise2.align.localmd(
                        refItem[1],
                        str(each.seq),
                        3, -1, -2, -1, -6, -4)
                aln = sorted(aln, key=operator.itemgetter(2))[0]
                print(refItem[0] + ' : ' + str(aln.score))
                if aln.score > maxScore:
                    maxScore = aln.score
                    finalAln = aln
                    gen = refItem[0]
            outDict.update({
                    each.id : {
                        'aln' : finalAln,
                        'genotype' : gen
                        }})
        print(outDict)


if __name__ == '__main__':
    refs = hbvRefs('references.fa')
    refs.throwOn('sanger.fa')
