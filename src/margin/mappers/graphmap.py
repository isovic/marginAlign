from margin.mappers.abstractMapper import AbstractMapper
from sonLib.bioio import system
import os

class GraphMap(AbstractMapper):
    def run(self, args=" -C "):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because GraphMap builds these crufty index files, copy to a temporary directory
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile))
        system("graphmap -I -r %s" % localReferenceFastaFile)
        system("graphmap %s -r %s -d %s -o %s" % (args, localReferenceFastaFile, self.readFastqFile, self.outputSamFile))

class GraphMapChain(GraphMap):
    def run(self):
        GraphMap.run(self)
        self.chainSamFile()

class GraphMapRealign(GraphMap):
    def run(self):
        GraphMap.run(self)
        self.realignSamFile()

class GraphMapAnchor(AbstractMapper):
    def run(self, args=" -C -a anchor "):
        localReferenceFastaFile = os.path.join(self.getLocalTempDir(), "ref.fa") #Because GraphMap builds these crufty index files, copy to a temporary directory
        system("cp %s %s" % (self.referenceFastaFile, localReferenceFastaFile))
        system("graphmap -I -r %s" % localReferenceFastaFile)
        system("graphmap %s -r %s -d %s -o %s" % (args, localReferenceFastaFile, self.readFastqFile, self.outputSamFile))

class GraphMapAnchorChain(GraphMapAnchor):
    def run(self):
        GraphMapAnchor.run(self)
        self.chainSamFile()

class GraphMapAnchorRealign(GraphMapAnchor):
    def run(self):
        GraphMapAnchor.run(self)
        self.realignSamFile()
