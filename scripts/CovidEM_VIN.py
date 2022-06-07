import numpy as np
from numpy.random import default_rng
import itertools
import logging

from collections import defaultdict
from collections import Counter

import sys
import argparse
import re 

#get special functions
from scipy.special import psi as digamma
from scipy.special import gammaln as lngamma
from itertools import compress

def expNormLogProb(logProbs):

    maxP = np.max(logProbs)
    
    ds = logProbs - maxP

    probs = np.exp(ds)
    
    probs /= probs.sum()
    
    return probs

MIN_ABUND = 0.025

MIN_AMP   = 0.5


class MixtureEM_VI():
    """ Class for mixture model over references using variational inference"""    
 
    DEFAULT_INIT_DELTA = 0.01
    DEFAULT_PI_PRIOR   = 0.001
    
    MIN_Z = 1.0e-9
 
    def __init__(self, N, R, A, mismatchMatrix, matchMatrix, mapA, piPrior = DEFAULT_PI_PRIOR, initDelta = DEFAULT_INIT_DELTA):
    
        self.R = R  # No. of references
    
        self.N = N #No. of reads
    
        self.A = A#No. of amplicons
    
        self.mismatchMatrix = mismatchMatrix
        
        self.matchMatrix = matchMatrix
        
        self.mapA = mapA # one hot mapping N -> A
        
        self.initDelta = initDelta
        
        self.deltaPC = 100.
        
        self.deltaPrior = self.initDelta*self.deltaPC
        
        self.deltaOnePrior = (1.0 - self.initDelta)*self.deltaPC
    
        self.piPrior = piPrior
    
        self.logExpDelta    = np.zeros(self.A)
        
        self.logExpDelta.fill(np.log(self.initDelta))
        
        self.logExpOneDelta = np.zeros(self.A)
         
        self.logExpOneDelta.fill(np.log(1.0 - self.initDelta))
    
        self.logExpPi = np.zeros(R)
        
        self.alphaPi = np.zeros(R)
        
        self.expPi = np.zeros(R)

        self.logExpPi.fill(np.log(self.piPrior/R))
    
        self.expZ = np.zeros((N,R))
    
    
    def update(self, maxIter = 100, minChange = 1.0e-3):
    
        nIter = 0
        lastElbo = float('inf') 
        fChange = float('inf') 
        while nIter < maxIter and fChange > minChange:
    
            #Update Z
            self.updateZ()
            
            self.updatePi()
            
            self.updateDelta()
            #Update delta
            
            elbo = self.calcELBO()
            
            fChange = abs(elbo - lastElbo)
            
            fString = str(nIter) + ',' + str(elbo) + ',' + str(fChange)
            
            logging.info(fString)
    
            lastElbo = elbo
            nIter += 1


    def updateZ(self):
    
        self.logRho = self.mismatchMatrix*np.dot(self.mapA,self.logExpDelta)[:,np.newaxis]
        
        self.logRho += self.matchMatrix*np.dot(self.mapA,self.logExpOneDelta)[:,np.newaxis]
        
        self.logRho += self.logExpPi[np.newaxis,:] 
    
        maxP = np.max(self.logRho,axis=1)
    
        ds = self.logRho - maxP[:,np.newaxis]

        expZ = np.exp(ds)
    
        row_sums = expZ.sum(axis=1)

        self.expZ = expZ / row_sums[:, np.newaxis]
        
    def updatePi(self):
    
        #Update Pi
        
        self.alphaPi = self.piPrior + np.sum(self.expZ,axis=0)        
        
        self.alphaSum = np.sum(self.alphaPi)
        
        self.logExpPi = digamma(self.alphaPi) - digamma(self.alphaSum)
        
        self.expPi = self.alphaPi/self.alphaSum
    
    def getVarPi(self):    
    
        varPi = self.expPi*(1.0 - self.expPi)
        
        varPi /= 1.0 + self.alphaSum
        
        return varPi
    
    def updateDelta(self):
        
        tM = self.expZ*self.matchMatrix
        
        self.mMatch = np.sum(np.dot(np.transpose(tM),self.mapA),axis=0) + self.deltaOnePrior
       
        tM = self.expZ*self.mismatchMatrix
        
        self.mMisMatch = np.sum(np.dot(np.transpose(tM),self.mapA),axis=0) + self.deltaPrior
    
        self.mTotal = self.mMatch + self.mMisMatch
    
        self.expDelta = self.mMisMatch/self.mTotal
    
        self.logExpDelta    = digamma(self.mMisMatch) - digamma(self.mTotal)
        
        self.logExpOneDelta =  digamma(self.mMatch) - digamma(self.mTotal)
    
    
    def logB(self,alpha):
    
        alphaSum = np.sum(alpha)
        
        return lngamma(alphaSum) - np.sum(lngamma(alpha))
    
    def calcELBO(self):
    
        tM = np.dot(self.mapA, self.logExpDelta)
    
        expLogLike = np.sum(self.expZ*self.mismatchMatrix*tM[:,np.newaxis]) 
    
        tM = np.dot(self.mapA, self.logExpOneDelta)
        
        expLogLike += np.sum(self.expZ*self.matchMatrix*tM[:,np.newaxis])
        
        expLogZ = np.sum(np.dot(self.expZ,self.logExpPi))
        
        expLogPi = lngamma(self.piPrior*self.R) - self.R*lngamma(self.piPrior) + (self.piPrior - 1.)*np.sum(self.logExpPi)
        
        expLogDelta = self.A*(lngamma(self.deltaPC) - lngamma(self.deltaPrior) - lngamma(self.deltaOnePrior)) + np.sum(self.deltaPrior*self.logExpDelta) + np.sum(self.deltaOnePrior*self.logExpOneDelta)
        
        minZ = self.expZ[self.expZ > MixtureEM_VI.MIN_Z]
        
        expLogQZ = np.sum(minZ*np.log(minZ))
        
        expLogQPi = np.sum((self.alphaPi - 1)*self.logExpPi) + self.logB(self.alphaPi)
        
        expLogQDelta = np.sum((self.mMatch - 1.0)*self.logExpOneDelta) +  np.sum((self.mMisMatch - 1.0)*self.logExpDelta) 
        
        expLogQDelta += np.sum(lngamma(self.mMatch + self.mMisMatch)) - np.sum(lngamma(self.mMatch) + lngamma(self.mMisMatch))
        
        #expLogDelta = lngamma(self.piPrior*self.R) - self.R*lngamma(self.piPrior) + (self.piPrior - 1.)*np.sum(self.logExpPi)
        
        ELBO = expLogLike + expLogZ + expLogPi + expLogDelta - expLogQZ - expLogQPi - expLogQDelta
    
        return ELBO
    
    
    def getAmpliconWeights(self):
    
        ampPi = np.dot(np.transpose(self.expZ),self.mapA)
    
        return ampPi
    
    def writeOutput(self,output_stub,ampl_list,ref_list,read_list):
    
        ampPi = self.getAmpliconWeights()
        
        aString = ','.join(ampl_list)
    
        aOutFile = output_stub + "_amp_pi_est.csv" 

        with open(aOutFile,'w') as f:
            print('Idx,Ref,' + aString,file=f)
        
            for r in range(self.R):
                zList = [str(x) for x in ampPi[r,:].tolist()]
                zString = ','.join(zList)        

                print('%d,%s,%s' % (r,ref_list[r],zString),file=f)
       

        rOutFile = output_stub + "_pi_est.csv" 

        sdPi = np.sqrt(self.getVarPi())
        with open(rOutFile,'w') as f:
            print('Idx,Ref,MeanFreq,StdFreq',file=f)
            for r in range(self.R):
                print('%d,%s,%f,%f' % (r,ref_list[r],self.expPi[r],sdPi[r]),file=f)
            
   
        dOutFile = output_stub + "_delta_est.csv"
   
        with open(dOutFile,'w') as f:
            print('Amp,MeanDelta',file=f)
            for n in range(self.A):
                print('%d,%s,%f' % (n,ampl_list[n],self.expDelta[n]),file=f)
   
        
        zOutFile = output_stub + "_z_est.csv"
   
        rString = ','.join(ref_list)
    
        with open(zOutFile,'w') as f:
            print('N,R,' + rString,file=f)
            for n in range(self.N):
                zList = [str(x) for x in self.expZ[n,:].tolist()]
                zString = ','.join(zList)
                print('%d,%s,%s' % (n,read_list[n],zString),file=f)
            
#MIN_AMP_FREQ = 10
    
def main(argv):
    parser = argparse.ArgumentParser()

    parser.add_argument("blast_file_m", help="merged alignments in tab format")
    
    parser.add_argument("read_lengths_m", help="read lengths tab delimited")
    
    parser.add_argument("mapping_file", help="mapping of amplicons to genomes")
    
    parser.add_argument("output_stub", help="file stub for output")

    parser.add_argument("-v","--verbose",action='store_true',help='Write logging info to stdout and log file')
    
    parser.add_argument("-r","--ref_file", default=None,help='File with list of reference sequences to use')
    
    parser.add_argument("-a","--amplicon_file", default=None,help='File with list of amplicon variants to use')
    
    parser.add_argument("-m","--max_iter", default=500,type=int,help = 'No. of iterations')
    
    parser.add_argument("-i","--min_read_length", default=0,type=int,help = 'Minimum merged read length for alignment to be used')
    
    args = parser.parse_args()
    
    #Set up logger
   # import ipdb; ipdb.set_trace()
    logging.basicConfig(filename=args.output_stub+'.log',  filemode='w',
                            format='%(asctime)s,%(msecs)d %(name)s %(levelname)s %(message)s',
                            datefmt='%H:%M:%S',
                            level=logging.DEBUG)
    
    
    if args.verbose:
        logging.getLogger().addHandler(logging.StreamHandler())
    
    var_genomes = defaultdict(list)
    
    
    ref_filter = set()
    bRefFilter = False
    if args.ref_file is not None:
        logging.info('Read reference filter file')
        bRefFilter = True
        with open(args.ref_file,'r') as f:
        
            for line in f:
        
                line = line.rstrip()
    
                ref_filter.add(line)
    else:
        logging.info('Not using reference filter file')
        
    bAmpliconFilter = False
    amplicon_filter = set()
    if args.amplicon_file is not None:
        logging.info('Read amplicon filter file')
        bAmpliconFilter = True
        with open(args.amplicon_file,'r') as f:
        
            for line in f:
        
                line = line.rstrip()
    
                amplicon_filter.add(line)
    else:
        logging.info('Not using amplicon filter file')
    
    var_filter = set()
    genome_amp_map = defaultdict(dict)
    logging.info('Reading amplicon mapping file')
    with open(args.mapping_file,'r') as f:
        
        for line in f:
        
            line = line.rstrip()
    
        
            toks=line.split('\t')
            
            var_amp = toks[1].split('_')[1]
            
            if toks[2] in ref_filter or not bRefFilter:
                var_genomes[toks[1]].append(toks[2])
                genome_amp_map[var_amp][toks[2]] = toks[1]
                if var_amp in amplicon_filter or not bAmpliconFilter:
            
                    var_filter.add(toks[1])
    
    
    read_lengths = {}
    
    logging.info('Reading read length file')
    with open(args.read_lengths_m,'r') as f:
        
        for line in f:
        
            line = line.rstrip()
    
        
            toks=line.split('\t')
    
            read_lengths[toks[0]] = int(toks[1])
    

    queries = set()
    
    errorDict = defaultdict(dict) 
    
    ampReadMap = {}
    
    readMaps = defaultdict(dict)
    
    refs = set()
    subjAmps = set()
    ampFreq = Counter()
    b = 0
    logging.info('Reading blast alignment file')
    with open(args.blast_file_m,'r') as f:
        
        for line in f:
        
            line = line.rstrip()
        
            if not line.startswith('#'):
        
                toks=line.split('\t')
            
                query = toks[0]
            
                subj = toks[1]
             
                length = int(toks[12])

                pid = 0.01*float(toks[3])
                
                mismatch = int(round(length*(1 - pid)))# int(toks[4]) + int(toks[5])
            
                subjectAmp = subj.split('_')[1]
            
                if length > args.min_read_length: 
                    
                    if subjectAmp in amplicon_filter or not bAmpliconFilter:
                        
                        if subj in var_filter or not bRefFilter:
                            ref = var_genomes[subj]
                        
                            refs.update(ref)
                            
                            if query not in queries:
                                    
                                queries.add(query)
                                
                                ampReadMap[query] = subjectAmp
                                
                                ampFreq[subjectAmp] += 1   
                                
                                subjAmps.add(subjectAmp)
                                    
                            readMaps[query][subj] = (mismatch,length - mismatch)
                    
                            queryAlign = toks[23]
                            subjectAlign = toks[24]
                    
                            
                    
            if (b % 10000 == 0):
                logging.info('Read ' + str(b) + ' alignments')
            
            b = b + 1
    
    #ampFilter = set( [k for k, v in ampFreq.items() if v > MIN_AMP_FREQ])
    ampFilter = subjAmps
    genomeMap = defaultdict(str)
    
    for g in refs:
        for amp in ampFilter:
            if g in genome_amp_map[amp]:
                genomeMap[g] += genome_amp_map[amp][g] + '$'
            else:
                genomeMap[g] += '*$'
    
    
    uniqueKeys = set(genomeMap.values())
    
    uMap = defaultdict(list)
    
    for (g,gkey) in genomeMap.items():
        uMap[gkey].append(g) 
    
    gChar = {}
    mapBack = {}
    for (mapg,gSet) in uMap.items():
        
        gSet = sorted(gSet)
    
        gChar[mapg] = gSet[0]
        mapBack[gSet[0]] = uMap[mapg]
    
    gChar = sorted(gChar.values())
    
    ref_list = sorted(gChar)
    
    R = len(ref_list) # number of references
    
    ref_map = {ref:i for (i,ref) in enumerate(ref_list)}
    
    mapFile = args.output_stub + "ref_map.csv"
    
    with open(mapFile,'w') as f:
        print('n,ID,NG,Genomes',file=f)
       
        for r in range(R):
            mString = '$'.join(mapBack[ref_list[r]])
            print('%d,%s,%d,%s' % (r,ref_list[r],len(mapBack[ref_list[r]]),mString),file=f)
    
    
    read_list = sorted(list(queries))
    
    N = len(queries) # number of reads
    print("queries %s"%N)
    ampl_list = sorted(list(subjAmps))
    
    A = len(ampl_list)
    
    ampl_map = {ref:i for (i,ref) in enumerate(ampl_list)}
    
    ampMap = np.zeros((N,A),dtype=int)
    
    for n in range(N):
        aname = ampReadMap[read_list[n]]
        
        ampMap[n,ampl_map[aname]] = 1
                        
    matchMatrix = np.zeros((N,R),dtype=int) #default no matches

    mismatchMatrix = np.zeros((N,R),dtype=int) #default all mismatches
    
    hitSlices = {}
    
    for r in range(N):
        
        read_id = read_list[r]
        
        rlen = read_lengths[read_id] 
        
        mismatchMatrix[r,:] = rlen
        
        hitSlice = []
        for hit,  (mismatch, mmatch) in readMaps[read_id].items():
        
            hGenomes = var_genomes[hit]
            
            for g in hGenomes: 
                if g in ref_list:
                    gidx = ref_map[g]
            
                    mismatchMatrix[r,gidx] = mismatch
                    matchMatrix[r,gidx] = mmatch
                    hitSlice.append(gidx)
            
        hitSlices[r] = tuple(sorted(hitSlice))
    
        if (r % 10000 == 0):
            logging.info('Processed ' + str(r) + ' reads')
    
    
    
    mixtureEM =  MixtureEM_VI(N, R, A, mismatchMatrix, matchMatrix, ampMap)
    
    logging.info('Running variational inference')
    mixtureEM.update(args.max_iter)
    
    logging.info('Writing output from VI')
    mixtureEM.writeOutput(args.output_stub,ampl_list,ref_list,read_list)
   
    aOutFile = args.output_stub + "amp_map.csv"
    
    with open(aOutFile,'w') as f:
        print('n,ID,Amp',file=f)
        for n in range(N):
            print('%d,%s,%s' % (n,read_list[n],ampReadMap[read_list[n]]),file=f)
   
    #Now rerun post filtering 
   
    logging.info('Filtering with min abundance: ' + str(MIN_ABUND))
    FilterRefs = mixtureEM.expPi > MIN_ABUND
   
    ampPi = mixtureEM.getAmpliconWeights()
   
    fAmpPi  = ampPi/ampPi.sum(axis=0,keepdims=1)
   
    FilterA = np.sum(fAmpPi > MIN_ABUND,axis=1) > MIN_AMP*A
   
    FilterRefs = np.logical_and(FilterRefs,FilterA)
    
    mismatchMatrix_Filt = mismatchMatrix[:,FilterRefs] 
    
    matchMatrix_Filt = matchMatrix[:,FilterRefs] 
    
    F = np.sum(FilterRefs)
    logging.info(str(F) + ' genomes after filtering')
    
    if F > 0:
        
    
        filt_ref_list = list(compress(ref_list, FilterRefs.tolist()))
    
        mixtureEM_Filt =  MixtureEM_VI(N, F, A, mismatchMatrix_Filt, matchMatrix_Filt, ampMap)
    
        logging.info('Run VI post filtering')
        
        mixtureEM_Filt.update(args.max_iter)
    
        mixtureEM_Filt.writeOutput(args.output_stub + '_Filt',ampl_list,filt_ref_list,read_list)
    
    else:
        logging.info('No genomes above minimum abundance not running filtering')
    
    
    logging.info('Program completed')


if __name__ == "__main__":
    main(sys.argv[1:])