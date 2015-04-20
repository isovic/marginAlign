The marginAlign package can be used to align reads to a reference genome and call single nucleotide variations (SNVs). It is specifically tailored for Oxford Nanopore Reads.

The package comes with three programs: marginAlign, a read aligner, marginCaller, a program to call single nucleotide variations, and marginStats, a program to compute simple qc stats on a sam file (alignment identity, coverage, insertion and deletion rates).

### Requirements
* git
* python 2.7
* pip/virtualenv (see below)

### Installation
To install the code run:

    git clone git://github.com/benedictpaten/marginAlign.git
    cd marginAlign
    git pull
    git submodule update --init
    make

This will build the code. Three executables: "marginAlign", "marginCaller" and "marginStats" are in the base directory
of the package. Place these binaries on your path if you wish to use them without referring to their absolute
path on the filesystem.

### Creating a virtual environment to handle python dependencies
marginAlign uses numpy, pysam and PyVcf (see requirements.txt for versions used). The system python could be used if these dependencies are present. Otherwise, a virtual environment can be created by running the following command in the margin base directory:

    virtualenv --no-site-packages --distribute env && source env/bin/activate && pip install -r requirements.txt

Creation of a virtual environment requires that machine has an existing installation of [pip](https://pip.pypa.io/en/latest/index.html) and [virtualenv](https://virtualenv.pypa.io/en/latest/index.html).

### Testing
To test the installation run:

    make test
    
This will run demo sequences through marginAlign and marginVar.
    
### Updating the installation
To update a marginAlign installation, from the base directory type:

    git pull
    git submodule update --init
    make clean
    make
    #If you're using virtualenv run the following command to update the python dependencies:
    virtualenv --no-site-packages --distribute env && source env/bin/activate && pip install -r requirements.txt
    

### jobTree

Both marginAlign and marginVar are [jobTree](https://github.com/benedictpaten/jobTree) scripts. This allows you to parallelise the execution of the scripts either using a single multi-processor machine or a cluster. The most important thing to note is that a jobTree script creates a "jobTree", a directory of files that manages the running of the job. As a result both scripts take a --jobTree DIRECTORY argument (by default "jobTree" in the current working directory if not specified). This directory should not exist prior to executing the script. It will be not be deleted after the script has finished (as it can be used to interogate the details of a run), and should be manually deleted.

### Running marginAlign

To access the help message type:

    marginAlign --help

or just:

    marginAlign

This will give you a full list of options. Most are related to the [jobTree](https://github.com/benedictpaten/jobTree) options, which control how the script is executed. 

To align a FASTQ file ("input.fastq") to a reference fasta file ("reference.fasta") and output the alignment in SAM format with marginAlign ("output.sam") using the "./jobTree" jobTree to manage the batch do:

    marginAlign input.fastq reference.fasta output.sam --jobTree ./jobTree

After executing the "./jobTree" directory should be deleted. This directory contains details of the run - it must be deleted before starting a new run that uses the same jobTree directory.

To enable EM training do, putting the trained model file in "output.hmm" do:

    marginAlign input.fastq reference.fasta output.sam --em --outputModel output.hmm --jobTree ./jobTree

To use a different model to the default one do:

    marginAlign input.fastq reference.fasta output.sam --inputModel input.hmm --jobTree ./jobTree

marginAlign expects the headers in the FASTQ file to be unique. Alternatively, a FASTQ file with unique headers can be created using the uniquifyFastq utility. To use uniquifyFastq do:

    uniquifyFastq input.fastq input_with_unique_headers.fastq

### Running marginCaller

To call single nucleotide variations from an existing alignment ("input.sam") sam file and the output ("output.vcf"):

    marginCaller input.sam reference.fa output.vcf --jobTree ./jobTree

To NOT marginalise over the read alignments do (this will just use the existing alignment, and will be much quicker):

    marginCaller input.sam reference.fa output.vcf --noMargin --jobTree ./jobTree

### Running marginStats

To calculate the median/avg/min/max identity of reads in a sam file do:
    
    marginStats input.sam read.fastq reference.fasta --identity

Other flags (see help) can be used to calculate other stats.

### Citing marginAlign/marginCaller/marginStats

Please cite margin align as:

    "Nat Methods. 2015 Apr;12(4):351-6. doi: 10.1038/nmeth.3290. Epub 2015 Feb 16.
    Improved data analysis for the MinION nanopore sequencer.
    Jain M1, Fiddes IT1, Miga KH1, Olsen HE1, Paten B1, Akeson M1."

