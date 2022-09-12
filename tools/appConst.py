import pathlib 
from string import Template


toolsPath = str(pathlib.Path(__file__).parent.resolve())
rootPath = '/'.join(toolsPath.split('/')[:-1])
blastDBpath = rootPath + '/BLASTdb/hbv' 
blastOutCols = 'qseqid sseqid qstart qend sstart send qseq sseq evalue length mismatch positive gaps'
blastCmd = f'{toolsPath}/blastn'

argsCmd = {
        '-db' : blastDBpath,
        '-query' : Template('$filePath'),
        '-outfmt' : "6 " + blastOutCols,
        '-max_target_seqs' : '1'
            }

resultDict = {
            
    'result' : {
        
        'conclusion' : None,

    
        }


        }

