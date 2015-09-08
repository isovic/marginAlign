import pysam, sys, os
from jobTree.src.bioio import reverseComplement, fastaRead, system, fastaWrite, \
cigarRead, logger, nameValue
from margin.utils import *
from cPecan import cPecanEm
from cPecan.cPecanEm import Hmm, SYMBOL_NUMBER
import numpy as np

def mergeChainedAlignedSegments(chainedAlignedSegments, refSequence, readSequence):
    """Makes a single alignment for the given chained reads. Will soft soft clip
    the unaligned prefix and suffix of the readSequence.
    
    From doc on building pysam line
    a = pysam.AlignedSegment()
    a.query_name = "read_28833_29006_6945"
    a.query_sequence="AGCTTAGCTAGCTACCTATATCTTGGTCTTGGCCG"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 32
    a.mapq = 20
    a.cigar = ( (0,10), (2,1), (0,25) )
    a.mrnm = 0
    a.mpos=199
    a.isize=167
    a.qual="<<<<<<<<<<<<<<<<<<<<<:<9/,&,22;;<<<"
    a.tags = ( ("NM", 1),
               ("RG", "L1") )
    """
    cAR = pysam.AlignedSegment()
    aR = chainedAlignedSegments[0]
    cAR.query_name = aR.query_name
    
    #Parameters we don't and therefore set properly
    #cAR.flag = aR.flag
    #cAR.mapq = aR.mapq
    #cAR.mrnm = 0
    #cAR.mpos=0
    #cAR.isize=0
    #cAR.qual = "<" * len(readSequence)
    #cAR.tags = aR.tags 
    cAR.next_reference_id = -1
    cAR.reference_start = aR.reference_start #Reference start
    cAR.is_reverse = aR.is_reverse
    cAR.query_sequence = reverseComplement(readSequence) if cAR.is_reverse else readSequence
    cAR.reference_id = aR.reference_id
    cigarList = []
    pPos = aR.reference_start
    #Iterate from the other end of the sequence if reversed
    pQPos = -(len(readSequence)-1) if cAR.is_reverse else 0 
        
    for aR in chainedAlignedSegments:
        assert cAR.is_reverse == aR.is_reverse
        #Add a deletion representing the preceding unaligned reference positions
        assert aR.reference_start >= pPos
        if aR.reference_start > pPos:
            cigarList.append((2, aR.reference_start - pPos))
            pPos = aR.reference_start 
    
        #Add an insertion representing the preceding unaligned read positions
        #make it a soft clip if it is the first chained alignment
        qPos = getFirstNonClippedPositionInRead(aR, readSequence)
        assert qPos >= pQPos
        if qPos > pQPos:
            cigarList.append((4 if aR == chainedAlignedSegments[0] else 1, qPos - pQPos)) 
            pQPos = qPos
        
        #Add the operations of the cigar, filtering hard and soft clipping
        for op, length in aR.cigar:
            assert op in (0, 1, 2, 4, 5)
            if op in (0, 1, 2):
                cigarList.append((op, length))
            if op in (0, 2): #Is match or deletion
                pPos += length
            if op in (0, 1): #Is match or insertion
                pQPos += length
        
    assert pPos <= len(refSequence)
    
    #Set reference end coordinate (which is exclusive)
    #cAR.reference_end = pPos #We don't do this because it is set by cigar string
    
    #Now add any trailing, necessary soft clipping
    if cAR.is_reverse:
        assert pQPos <= 1
        if pQPos < 1:
            cigarList.append((4, -pQPos + 1))
    else:
        assert pQPos <= len(readSequence)
        if pQPos < len(readSequence):
            cigarList.append((4, len(readSequence) - pQPos))
    
    cAR.cigar = tuple(cigarList)
    
    #Check ops
    for op, length in cAR.cigar: #We should have no hard clipped ops
        assert op in (0, 1, 2, 4)
     
    #Reference sequence check coordinates
    assert sum([ length for op, length in cigarList if op in (0, 2)]) == cAR.reference_end - cAR.reference_start
    assert cAR.reference_start >= 0 and cAR.reference_start < len(refSequence)
    assert cAR.reference_end >= 0 and cAR.reference_end <= len(refSequence)
    
    #Read sequence check coordinates
    assert cAR.query_alignment_start >= 0 and cAR.query_alignment_start < len(readSequence)
    assert cAR.query_alignment_end >= 0 and cAR.query_alignment_end <= len(readSequence)
    assert cAR.query_alignment_start + sum([ length for op, length in cigarList if op in (0, 1)]) == cAR.query_alignment_end
    
    return cAR

def chainFn(alignedSegments, refSeq, readSeq, scoreFn=\
            lambda alignedSegment, refSeq, readSeq : \
            sum([ length for op, length in alignedSegment.cigar if op == 0]), maxGap=200):
     #Score function is number of aligned pairs
    """Gets the highest scoring chain of alignments on either the forward or reverse 
    strand. Score is (by default) number of aligned positions.
    """
    def getStartAndEndCoordinates(alignedSegment):
        """Gets the start and end coordinates in both the reference and query, using coordinates
        relative to the original read and reference equence
        """
        return alignedSegment.reference_start, getFirstNonClippedPositionInRead(alignedSegment, readSeq), \
        alignedSegment.reference_end-1, getLastNonClippedPositionInRead(alignedSegment, readSeq) 
    
    alignedSegmentToScores = dict([ (aR, scoreFn(aR, refSeq, readSeq)) for aR in alignedSegments])
    alignedSegmentToCoordinates = dict([ (aR, getStartAndEndCoordinates(aR)) for \
                                     aR in alignedSegments])
    alignedSegmentPointers = {}
    
    #Currently uses sloppy quadratic algorithm to find highest chain
    alignedSegments = sorted(alignedSegments, key=lambda aR : alignedSegmentToCoordinates[aR][0]) 
    #Sort by reference coordinate
    for i in xrange(len(alignedSegments)):
        aR = alignedSegments[i]
        rStart, qStart, rEnd, qEnd = alignedSegmentToCoordinates[aR]
        score = alignedSegmentToScores[aR]
        for j in xrange(i): #Look at earlier alignments in list
            aR2 = alignedSegments[j]
            rStart2, qStart2, rEnd2, qEnd2 = alignedSegmentToCoordinates[aR2]
            assert rStart2 <= rStart
            if rStart > rEnd2 and qStart > qEnd2 and aR.is_reverse == aR2.is_reverse and \
            rStart - rEnd2 + qStart - qEnd2 <= maxGap and \
            score + alignedSegmentToScores[aR2] > alignedSegmentToScores[aR]: 
            #Conditions for a chain
                alignedSegmentToScores[aR] = score + alignedSegmentToScores[aR2]
                alignedSegmentPointers[aR] = aR2
    
    #Now find highest scoring alignment 
    aR = max(alignedSegments, key=lambda aR : alignedSegmentToScores[aR])
    
    #Construct chain of alignedSegments
    chain = [ aR ]
    while aR in alignedSegmentPointers:
        aR = alignedSegmentPointers[aR]
        chain.append(aR)
    chain.reverse()
    
    return chain

def chainSamFile(samFile, outputSamFile, readFastqFile, referenceFastaFile, 
                 chainFn=chainFn):
    """Chains together the reads in the SAM file so that each read is covered by a 
    single maximal alignment
    """
    sam = pysam.Samfile(samFile, "r" )
    refSequences = getFastaDictionary(referenceFastaFile) #Hash of names to sequences
    
    alignmentsHash = {}
    for aR in samIterator(sam): #Iterate on the sam lines and put into buckets by read
        #This should be improved, because the whole sam file is being stored in memory
        if aR.query_name not in alignmentsHash:
            alignmentsHash[aR.query_name] = {}
        if aR.reference_id not in alignmentsHash[aR.query_name]:
            alignmentsHash[aR.query_name][aR.reference_id] = []
        alignmentsHash[aR.query_name][aR.reference_id].append(aR)

    #Now write out the sam file
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    
    #Chain together the reads
    chainedAlignedSegments = []
    for readName, readSeq, qualValues in fastqRead(readFastqFile):
        readName = readName.split()[0] #Remove any white space from the name
        if readName in alignmentsHash:
            for refID in alignmentsHash[readName].keys():
                alignedSegments = alignmentsHash[readName][refID]
                refSeq = refSequences[sam.getrname(refID)]
                chainedAlignedSegments.append(mergeChainedAlignedSegments(chainFn(alignedSegments, 
                                              refSeq, readSeq), refSeq, readSeq))
            alignmentsHash.pop(readName)
    assert len(alignmentsHash) == 0 #All reads in the sam file should be in the input sequence file
        
    #Sort chained alignments by reference coordinates
    chainedAlignedSegments.sort(key=lambda aR : (sam.getrname(aR.reference_id), \
                                                 aR.reference_start, aR.reference_end)) 
    
    for cAR in chainedAlignedSegments:
        outputSam.write(cAR)
    sam.close()
    outputSam.close()
    
def learnModelFromSamFileTargetFn(target, samFile, readFastqFile, 
                                  referenceFastaFile, options):
    """Does expectation maximisation on sam file to learn the hmm for the sam file.
    """
    #Get cigars and reads fasta file
    cigars = os.path.join(target.getGlobalTempDir(), "temp.cigar")
    fHCigars = open(cigars, 'w')
    reads = os.path.join(target.getGlobalTempDir(), "temp.fa")
    fHReads = open(reads, 'w')
    sam = pysam.Samfile(samFile, "r" )
    for aR, counter in zip(sam, xrange(sys.maxint)): #Iterate on the sam lines realigning them in parallel            
        aR.query_name = aR.query_name + "_%s" % counter
        fHCigars.write(getExonerateCigarFormatString(aR, sam) + "\n")
        fastaWrite(fHReads, aR.query_name, aR.seq)
    fHCigars.close(); fHReads.close()
    
    unnormalisedOutputModel = os.path.join(target.getGlobalTempDir(), 
                                           "unnormalisedOutputModel.hmm")
    target.addChildTargetFn(cPecanEm.expectationMaximisationTrials, 
                            args=(" ".join([reads, referenceFastaFile ]), cigars, 
                                  unnormalisedOutputModel, options))
    
    #Now set up normalisation
    target.setFollowOnTargetFn(learnModelFromSamFileTargetFn2, 
                               args=(unnormalisedOutputModel, options))

def learnModelFromSamFileTargetFn2(target, unnormalisedOutputModel, options):
    hmm = Hmm.loadHmm(unnormalisedOutputModel)
    setHmmIndelEmissionsToBeFlat(hmm)
    #Normalise background emission frequencies, if requested to GC% given
    normaliseHmmByReferenceGCContent(hmm, 0.5)
    hmm.write(options.outputModel)
    
toMatrix = lambda e : map(lambda i : e[SYMBOL_NUMBER*i:SYMBOL_NUMBER*(i+1)], 
                          xrange(SYMBOL_NUMBER))
fromMatrix = lambda e : reduce(lambda x, y : list(x) + list(y), e)
    
def normaliseHmmByReferenceGCContent(hmm, gcContent):
    """Normalise background emission frequencies to GC% given
    """
    for state in range(hmm.stateNumber):
        if state not in (2, 4): #Don't normalise GC content of insert states 
            #(as they don't have any ref bases!)
            n = toMatrix(hmm.emissions[(SYMBOL_NUMBER**2) * 
                                       state:(SYMBOL_NUMBER**2) * (state+1)])
            hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
            fromMatrix(map(lambda i : map(lambda j : (n[i][j]/sum(n[i])) * 
            (gcContent/2.0 if i in [1, 2] else (1.0-gcContent)/2.0), range(SYMBOL_NUMBER)), 
                           range(SYMBOL_NUMBER))) #Normalise

def setHmmIndelEmissionsToBeFlat(hmm):
    """Set indel emissions to all be flat
    """
    for state in range(1, hmm.stateNumber):
        hmm.emissions[(SYMBOL_NUMBER**2) * state:(SYMBOL_NUMBER**2) * (state+1)] = \
        [1.0/(SYMBOL_NUMBER**2)]*SYMBOL_NUMBER**2  

def modifyHmmEmissionsByExpectedVariationRate(hmm, substitutionRate):
    #Normalise background emission frequencies, if requested to GC% given
    n = toMatrix(map(lambda i : (1.0-substitutionRate) if i % SYMBOL_NUMBER == \
                     i / SYMBOL_NUMBER else substitutionRate/(SYMBOL_NUMBER-1), 
                     xrange(SYMBOL_NUMBER**2)))
    hmm.emissions[:SYMBOL_NUMBER**2] = fromMatrix(np.dot(toMatrix(hmm.emissions[:SYMBOL_NUMBER**2]), n))

def realignSamFileTargetFn(target, samFile, outputSamFile, readFastqFile, 
                           referenceFastaFile, options, chainFn=chainFn):
    """Chains and then realigns the resulting global alignments, using jobTree to 
    do it in parallel on a cluster.
    Optionally runs expectation maximisation.
    """
    #Optionally chain the sam file
    if not options.noChain:
        target.logToMaster("Going to chain sam file: %s" % samFile)
        tempSamFile = os.path.join(target.getGlobalTempDir(), "temp.sam")
        chainSamFile(samFile, tempSamFile, readFastqFile, referenceFastaFile, chainFn)
        samFile = tempSamFile
    
    #If we do expectation maximisation we split here:
    if options.em:
        target.logToMaster("Going to run EM training with sam file: %s, read fastq file: %s, \
        reference fasta file: %s" \
                           % (samFile, readFastqFile, referenceFastaFile))
        target.addChildTargetFn(learnModelFromSamFileTargetFn, args=(samFile, 
                                    readFastqFile, referenceFastaFile, options))

    options.hmmFile = options.outputModel if options.em else options.inputModel #This
    #setups the hmm to be used the realignment function
    
    target.logToMaster("Going to realign sam file: %s to create output sam file: %s \
    with match gamma %s and gap gamma %s and model %s" % (samFile, outputSamFile, 
                options.gapGamma, options.matchGamma, options.hmmFile))
    
    target.setFollowOnTargetFn(paralleliseSamProcessingTargetFn, 
                               args=(samFile, 
                                     referenceFastaFile, outputSamFile, 
                                     realignCigarTargetFn, realignSamFile3TargetFn,
                                     options))
    
def realignCigarTargetFn(target, exonerateCigarStringFile, referenceSequenceName, 
                         referenceSequence, querySequenceFile, 
                         outputCigarFile, options):
    #Temporary files
    tempRefFile = os.path.join(target.getLocalTempDir(), "ref.fa")
    tempReadFile = os.path.join(target.getLocalTempDir(), "read.fa")
    
    #Write the temporary reference file.
    fastaWrite(tempRefFile, referenceSequenceName, referenceSequence) 
    
    #For each cigar string
    for exonerateCigarString, (querySequenceName, querySequence) in \
    zip(open(exonerateCigarStringFile, "r"), fastaRead(querySequenceFile)):
        fastaWrite(tempReadFile, querySequenceName, querySequence)
        #Call to cPecanRealign
        loadHmm = nameValue("loadHmm", options.hmmFile)
#        system("echo %s | cPecanRealign %s %s --diagonalExpansion=10 \
 #       --splitMatrixBiggerThanThis=3000 %s --gapGamma=%s --matchGamma=%s >> %s" % \
  #             (exonerateCigarString[:-1], tempRefFile, tempReadFile, loadHmm, 
   #             options.gapGamma, options.matchGamma, outputCigarFile))
        try:
            command = "echo %s | cPecanRealign %s %s --diagonalExpansion=10 \
            --splitMatrixBiggerThanThis=3000 %s --gapGamma=%s --matchGamma=%s >> %s" % \
                   (exonerateCigarString[:-1], tempRefFile, tempReadFile, loadHmm, 
                    options.gapGamma, options.matchGamma, outputCigarFile);
            system(command)
            target.logToMaster('[good] ' + command + '\n');
        except Exception, e:            
            target.logToMaster('Caught an exception! qname = "%s"\n' % querySequenceName);
            target.logToMaster('len(exonerateCigarString[:-1]) = %d\n' % (len(exonerateCigarString[:-1])));
            target.logToMaster('[bad] Command that caused the exception:\n');
            target.logToMaster("echo %s | cPecanRealign %s %s --diagonalExpansion=10 --splitMatrixBiggerThanThis=3000 %s --gapGamma=%s --matchGamma=%s >> %s" % (exonerateCigarString[:-1], tempRefFile, tempReadFile, loadHmm, options.gapGamma, options.matchGamma, outputCigarFile));
            target.logToMaster('\n');
            target.logToMaster('\n');
            target.logToMaster(str(e) + '\n');
            target.logToMaster('\n');
            continue;

def realignSamFile3TargetFn(target, samFile, referenceFastaFile, 
                            outputSamFile, tempCigarFiles, options):
    #Setup input and output sam files
    sam = pysam.Samfile(samFile, "r" )
    
    #Replace the cigar lines with the realigned cigar lines
    outputSam = pysam.Samfile(outputSamFile, "wh", template=sam)
    def cigarIterator():
        #Iterates over all the cigars in the temp files.
        for tempCigarFile in tempCigarFiles:
            for pA in cigarRead(open(tempCigarFile)):
                yield pA 
        yield None #This is put in to cause an error if there is fewer 
        #cigars than pairwise alignments
    for aR, pA in zip(samIterator(sam), cigarIterator()): #Iterate on the sam lines 
        #realigning them in parallel
        
        #Replace alignment by converting exonerate ops to aligned read ops,
        #adding soft clipping/hard clipping of unaligned prefix and suffix of read
        ops = [ ]
        if len(aR.cigar) > 0 and aR.cigar[0][0] == 5:
            ##Add any hard clipped prefix
            ops.append(aR.cigar[0])
        if aR.query_alignment_start > 0:
            ops.append((4, aR.qstart))
        ops += map(lambda op : (op.type, op.length), pA.operationList)
        if aR.query_alignment_end < len(aR.query_sequence):
            ops.append((4, len(aR.query_sequence) - aR.query_alignment_end))
        if len(aR.cigar) > 1 and aR.cigar[-1][0] == 5: 
            ##Add any hard clipped suffix
            ops.append(aR.cigar[-1])
        
        #Checks the final operation list 
        ##Correct for the read
        assert sum(map(lambda (type, length) : length if type in (0,1,4) else 0, ops)) == \
        sum(map(lambda (type, length) : length if type in (0,1,4) else 0, aR.cigar))
        ##Correct for the reference
        assert sum(map(lambda (type, length) : length if type in (0, 2) else 0, ops)) == \
        aR.reference_end - aR.reference_start
        
        aR.cigar = tuple(ops)
        
        #Write out
        outputSam.write(aR)
    
    #Finish up
    sam.close()
    outputSam.close()
    
    target.logToMaster("Realigned sam file: %s to create output sam file: %s" % \
                       (samFile, outputSamFile))
