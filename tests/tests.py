import unittest
import time
import os, sys, numpy, pysam
from margin.marginCallerLib import vcfRead
from margin.utils import ReadAlignmentStats, pathToBaseNanoporeDir, getFastaDictionary
from cPecan.cPecanEm import Hmm
from sonLib.bioio import system, parseSuiteTestOptions, logger, getBasicOptionParser
import vcf
import numpy

"""Basic system level tests for marginAlign, marginCaller and marginStats scripts.
"""

longTests = False

class TestCase(unittest.TestCase):
    def setUp(self):
        self.marginAlign = self.getFile("marginAlign") #Path to marginAlign binary
        self.marginCaller = self.getFile("marginCaller") #Path to marginCall binary
        self.marginStats = self.getFile("marginStats") #Path to marginStats binary
        self.modifyHmm = self.getFile("scripts/modifyHmm")
        #The following are files used
        self.readFastqFile1 = self.getFile("tests/reads.fq") if longTests else \
        self.getFile("tests/lessReads.fq")
        self.referenceFastaFile1 = self.getFile("tests/references.fa")
        self.outputSamFile = self.getFile("tests/test.sam")
        self.inputHmmFile = self.getFile("tests/input.hmm")
        self.outputHmmFile = self.getFile("tests/output.hmm")
        self.inputSamFile1 = self.getFile("tests/input.sam")
        self.outputVcfFile = self.getFile("tests/output.vcf")
        self.jobTree = self.getFile("tests/testJobTree")
        
        #For testing margin caller
        self.mutationsFile = self.getFile("tests/mutations.txt")
        self.mutatedReferenceFastaFile = self.getFile("tests/referencesMutated.fa")
        self.inputSamFileForMutatedReferenceFile = self.getFile("tests/inputBigMutations.sam") #This is aligned against the mutated reference
        self.inputSamFileForMutatedReferenceFileLast = self.getFile("tests/inputBigMutationsLast.sam") 
        self.inputSamFileForMutatedReferenceFileBwa = self.getFile("tests/inputBigMutationsBwa.sam") 
        self.inputSamFileForMutatedReferenceFileGraphMap = self.getFile("tests/inputBigMutationsBwa.sam") 
        self.readFastqFile2 = self.getFile("tests/reads.fq")
        
        unittest.TestCase.setUp(self)

    def getFile(self, file):
        return os.path.join(pathToBaseNanoporeDir(), file)

    def tearDown(self):
        unittest.TestCase.tearDown(self)
        # Clean up
        system("rm -rf %s %s %s %s" % (self.outputSamFile, self.outputHmmFile, 
                                    self.outputVcfFile, self.jobTree))

    def validateSam(self, samFile, readFastqFile, referenceFastaFile):
        """Checks if a sam file is valid.
        """
        # Check if samfile exists
        self.assertTrue(os.path.isfile(samFile))
        #The call calculate identity will run a lot of internal consistency checks
        #as it calculates the alignment identity.
        return ReadAlignmentStats.getReadAlignmentStats(samFile, readFastqFile, 
                                                        referenceFastaFile, globalAlignment=True)
    
    def validateVcf(self, vcfFile, referenceFastaFile, mutationsFile):
        #Load reference sequences
        referenceSequences = getFastaDictionary(referenceFastaFile)
        #Load mutations
        mutations = set(map(lambda x : (x[0], int(x[1])+1, x[2]), \
                            map(lambda x : x.split(), open(mutationsFile, 'r'))))
        #Load VCF mutations
        imputedMutations = vcfRead(vcfFile)
            
        #print "Known mutations", sorted(list(mutations))
        #print "Imputed mutations", sorted(list(imputedMutations))
        
        #Compare mutation sets
        intersectionSize = float(len(mutations.intersection(imputedMutations)))
        #Return precision, recall, number of mutations called, number of known mutations 
        return intersectionSize/len(imputedMutations) if len(imputedMutations) else 0.0, \
    intersectionSize/len(mutations) if len(mutations) else 0.0, len(imputedMutations), len(mutations) 
    
    def checkHmm(self, hmmFile):
        Hmm.loadHmm(hmmFile) #This performs a bunch of internal consistency checks 
       
    def runMarginAlign(self, readFastqFile, referenceFastaFile, args=""):
        startTime = time.time()
        system("\t".join([ self.marginAlign, readFastqFile,
                         referenceFastaFile, self.outputSamFile, "--jobTree=%s" % self.jobTree, args ]))
        runTime = time.time() - startTime
        readAlignmentStats = self.validateSam(self.outputSamFile, readFastqFile, referenceFastaFile)
        #Get some stats to print
        identity = numpy.average(map(lambda rAS : rAS.identity(), readAlignmentStats))
        mismatchesPerAlignedBase = numpy.average(map(lambda rAS : rAS.mismatchesPerAlignedBase(), readAlignmentStats))
        insertionsPerReadBase = numpy.average(map(lambda rAS : rAS.insertionsPerReadBase(), readAlignmentStats))
        deletionsPerReadBase = numpy.average(map(lambda rAS : rAS.deletionsPerReadBase(), readAlignmentStats))
        
        logger.info("Ran marginAlign with args: %s, with reference: %s and reads: %s. \
        Got identity: %s, Mismatches per aligned base: %s, Insertions per read base: %s, \
        Deletions per read base: %s, Took: %s seconds" % \
                    (args, readFastqFile, referenceFastaFile, identity, 
                     mismatchesPerAlignedBase, insertionsPerReadBase,
                     deletionsPerReadBase, runTime))
        system("rm -rf %s" % self.jobTree)
    
    ###The following functions test marginAlign
    
    def testMarginAlignDefaults(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1)
    
    def testMarginAlignNoChain(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noChain")
    
    def testMarginAlignNoRealign(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noRealign")
    
    def testMarginAlignNoRealignNoChain(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noRealign --noChain")

    def testMarginAlignEm(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--em --outputModel %s" % self.outputHmmFile)
        self.checkHmm(self.outputHmmFile)
    
    def testMarginAlignEmNoChain(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--noChain --em --outputModel %s" % self.outputHmmFile)
        self.checkHmm(self.outputHmmFile)

    def testMarginAlignLoadCustomInputModel(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--inputModel %s" % self.inputHmmFile)

    def testMarginAlignBwa(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--bwa")
    
    def testMarginAlignBwaNoRealign(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--bwa --noRealign")
    
    def testMarginAlignGraphMap(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--graphmap")
    
    def testMarginAlignGraphMapNoRealign(self):
        self.runMarginAlign(self.readFastqFile1, self.referenceFastaFile1, "--graphmap --noRealign")
    
    #The following tests marginCaller
    
    def runMarginCaller(self, samFile, referenceFastaFile, mutationsFile, args=""):
        startTime = time.time()
        system("\t".join([ self.marginCaller, samFile, referenceFastaFile,
                         self.outputVcfFile, "--jobTree=%s" % self.jobTree, args ]))
        runTime = time.time() - startTime
        precision, recall, numberOfCalls, numberOfKnownMutations = self.validateVcf(self.outputVcfFile, 
                                             referenceFastaFile, mutationsFile)
        logger.info("Ran marginCaller with args: %s, with reference: %s and sam: %s. \
        Got: %s precision, Got: %s recall, Number of calls: %s, Number of known mutations: %s,\
        Took: %s seconds" % (args, samFile, referenceFastaFile, precision, recall, 
                             numberOfCalls, numberOfKnownMutations, runTime))
        system("rm -rf %s" % self.jobTree)
    
    def testMarginCallerDefaults(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFile, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile)
    
    def testMarginCallerNoMargin(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFile, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile, "--noMargin")
        
    def testMarginCallerNoMarginLast(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFileLast, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile)
    
    def testMarginCallerNoMarginBwa(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFileBwa, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile)

    def testMarginCallerNoMarginGraphMap(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFileGraphMap, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile)
        
    def testMarginCallerNoMarginLastNoMargin(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFileLast, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile, "--noMargin")
    
    def testMarginCallerNoMarginBwaNoMargin(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFileBwa, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile, "--noMargin")

    def testMarginCallerNoMarginGraphMapNoMargin(self):
        self.runMarginCaller(self.inputSamFileForMutatedReferenceFileGraphMap, 
                             self.mutatedReferenceFastaFile,
                             self.mutationsFile, "--noMargin")
    #Full integrative test that runs EM to train a model, then uses the resulting
    #model and alignment to calculate SNPs
    
    def testMarginAlignAndCallerTogether(self):
        if longTests:
            #This will run margin-align to get the output-model
            self.runMarginAlign(self.readFastqFile2, self.mutatedReferenceFastaFile, \
                                args="--em --outputModel %s" % self.outputHmmFile)
            #Establish the accuracy of the alignment
            self.runMarginCaller(self.outputSamFile, 
                                 self.mutatedReferenceFastaFile,
                                 self.mutationsFile, "--noMargin")
            #Establish the accuracy of the model before modification
            self.runMarginCaller(self.outputSamFile, 
                                 self.mutatedReferenceFastaFile,
                                 self.mutationsFile, "--alignmentModel=%s --errorModel=%s" % 
                                 (self.outputHmmFile, self.outputHmmFile))
            #Modify the HMM to make it more tolerant other substitutions / normalise the GC content / simplify emission probs
            system("%s %s %s --gcContent=0.5 --substitutionRate=0.2 --setFlatIndelEmissions" % (self.modifyHmm, self.outputHmmFile, self.outputHmmFile))
            #Establish the accuracy of the model after modification
            self.runMarginCaller(self.outputSamFile, 
                                 self.mutatedReferenceFastaFile,
                                 self.mutationsFile, "--alignmentModel=%s --errorModel=%s" % 
                                 (self.outputHmmFile, self.outputHmmFile))
    
    #This runs margin stats (just to ensure it runs without falling over)
    def testMarginStats(self):
        system("%s %s %s %s --identity --mismatchesPerAlignedBase --readCoverage \
        --deletionsPerReadBase --insertionsPerReadBase --printValuePerReadAlignment" % \
        (self.marginStats, self.inputSamFile1, self.readFastqFile1, self.referenceFastaFile1))
    
def main():
    parser = getBasicOptionParser()
    parser.add_option("--longTests", dest="longTests", action="store_true",
                      help="Run longer, more complete tests (with more reads)",
                      default=False)
    options, args = parseSuiteTestOptions(parser)
    global longTests
    longTests = options.longTests
    sys.argv = sys.argv[:1]
    unittest.main()
        
if __name__ == '__main__':
    main()