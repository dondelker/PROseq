#!/usr/bin/python3
# Writes out a bigWig file of coverage of a .bam file.
# This is along the lines of deepTools' bamCoverage, but
# with somewhat different options.

# Requires bedtools and samtools to be in $PATH.
# Assumes that the .bam file is already sorted.
# Also assumes that uniquely-mapped-read counts are available.

import argparse
import os
import os.path
import re
import subprocess
import sys
import tempfile

import pysam

import pdb

parser = argparse.ArgumentParser(description=
    'Writes coverage of a BAM file, in bigWig format.')
# For unstranded files, .bw will be appended.
# For stranded files, .plus.bw and .minus.bw will be appended.

parser.add_argument('BAM_input',
    help='BAM file to read')
parser.add_argument('bigWig_output',
    help='bigWig file to write out.')

parser.add_argument('--type',
    choices=['ChIP', 'PROseq', 'PROcap', 'RNA', '-5', '-3', 'paired'],
    default='ChIP', help=
'''type of coverage to compute. 'ChIP' will result in one
unstranded bigWig file, while 'PROseq', 'PROcap', 'RNA',
'-5', and '-3' will result in a pair of stranded bigWig files.
'-5' and '-3' will just get coverage at that end (thus '-3' is
a synonym for PROseq.)
'paired' will treat as stranded paired-end reads, and include
the fragment between the ends.''')

parser.add_argument('--unnormalized',
    default=False, action='store_true',
    help='if true, don\'t normalize to RPM')

args = parser.parse_args()

def getUniquelyAlignedReads(bamFile):
    """Gets count of unique aligned reads.

    The count is in a file that's generated when the reads are aligned.
    bamFile: name of the BAM file
    Returns: number of unique mapped reads (or None on error)."""
    with open(bamFile + '.stats', 'r') as f:
        # parse first line
        line = f.readline()
        m = re.match('^(?:Aligned reads: )?(\d+)$', line)
        if m:
            return int(m.group(1))
        if not re.match('^[Uu]nique.*reads:', line):
            return None
        # parse second line
        line = f.readline()
        m = re.match('^(\d+)$', line)
        if m:
            return int(m.group(1))
        return None

def countReads(bamFile):
    """Counts aligned reads, using samtools.

    This is intended to be the number of mapped reads (_not_ alignments).
    It seems to be at least very close to the number of read names in
    the .bam file.
    bamFile: name of the BAM file
    Returns: number of mapped reads (or None on error)."""
    output = subprocess.check_output(['samtools', 'view',
        # SAM/BAM spec. indicates there should only be one line with
        #   FLAG & 0x900 == 0
        # reads with 0x4 set are unmapped (although we're usually only
        # including mapped reads in BAM files anyway)
        '-F', str(0x904),
        '-c', bamFile])
    return int(output)

def writeSizes(bamFile):
    """Writes out the sizes of chromosomes, based on the BAM header.

    (Using pybedtools would be simpler in some ways).
    bamFile: name of the BAM file to get information for
    Returns: a file object, for a file of chromosome sizes
        (which I think will be deleted when the file object is
        garbage-collected)."""
    # get headers using samtools
    r = subprocess.run(['samtools', 'view', '-H', bamFile],
        stdout=subprocess.PIPE)
    # create temporary file
    f = tempfile.NamedTemporaryFile()
    # loop through lines, and get those including sequence info
    for line in r.stdout.decode().splitlines():
        # here, presumably the \d+ will greedily read the number,
        # whether or not there's anything after it
        m = re.match('^@SQ\tSN:(\S+)\tLN:(\d+)', line)
        if m:
            s = m.group(1) + '\t' + m.group(2) + '\n'
            f.write(s.encode())
    # flush the file, to make sure other programs can read it
    f.flush()
    return(f)

def writeBigWig(bamFile, bigWigFile, flags,
        normalize=True, strand=None, negate=False):
    """Writes coverage for one bigWig file.
    bamFile: the (sorted) BAM file to read
    bigWigFile: the bigWig file to write
    flags: flags to pass to 'bedtools genomecov'
    normalize: if True, normalize to RPM (reads per million)
    strand: if set to '+' or '-', only include reads
        from that strand
    negate: if True, then negate output
    Side effects: writes out a bigWig coverage file.
        The chromosome names and sizes are read from
        the BAM file.
        If the output file exists, doesn't do anything.
    """
    # if the input file isn't present, throw an error
    if not os.path.exists(bamFile):
        print('can\'t open file: ' + bamFile)
        sys.exit(1)
    # if the output file is present, don't do anything
    if os.path.exists(bigWigFile):
        return
    # make sure the directory exists
    outputDir = os.path.dirname(bigWigFile)
    if outputDir != '':
        os.makedirs(outputDir, exist_ok=True)
    # possibly compute normalization factor
    scale = 1
    if normalize:
        # for now, getting read counts using samtools
        readCount = countReads(bamFile)
        print('readcount = ' + str(readCount))
        scale = 1e6 / float(readCount)
    # possibly negate scale
    if negate:
        scale = -scale
    # get chromosome sizes files
    genomeFile = writeSizes(bamFile)
    # create a temporary bedGraph file
    tmpBedGraphFile = tempfile.NamedTemporaryFile(
        suffix='.bedGraph')
    # write out bedGraph file
    args = ['bedtools', 'genomecov', '-bga',
        '-ibam', bamFile,
        '-g', genomeFile.name,
        '-scale', str(scale)]
    if strand:
        args = args + ['-strand', strand]
    args = args + flags
    print(args)
    p1 = subprocess.call(args, stdout = tmpBedGraphFile)
    # also sort
    tmpBedGraphFile2 = tempfile.NamedTemporaryFile(
        suffix='.bedGraph')
    p2 = subprocess.call(['bedtools', 'sort', '-i', tmpBedGraphFile.name],
      stdout = tmpBedGraphFile2)
    # convert to bigWig format
    subprocess.call(['bedGraphToBigWig',
        tmpBedGraphFile2.name,
        genomeFile.name,
        bigWigFile])
    # ??? I think the temp file is removed automatically

# ChIP data: just write unstranded coverage
if args.type == 'ChIP':
    # jtb: now adding '-split' (to avoid issues from long
    # introns, although those should be avoided anyway)
    writeBigWig(args.BAM_input, args.bigWig_output + '.bw',
        ['-split'], normalize = not args.unnormalized)
    sys.exit(0)

# PROseq or RNAseq data: write stranded coverage in two
# separate files (the only difference is that for
# PROseq, use coverage of the 3' end, but for RNAseq,
# use exon coverage)
if args.type in ['PROseq', 'RNA', '-5', '-3']:
    # set the appropriate flags
    flags = ['-split']
    if args.type in ['-5']:
        flags = ['-5']
    if args.type in ['PROseq', '-3']:
        flags = ['-3']
    # write out the individual strands
    writeBigWig(args.BAM_input, args.bigWig_output + '.plus.bw',
        flags, normalize = not args.unnormalized,
        strand = '+', negate = False)
    # same, but negated
    writeBigWig(args.BAM_input, args.bigWig_output + '.minus.bw',
        flags, normalize = not args.unnormalized,
        strand = '-', negate = True)
    sys.exit(0)

# Lastly, and most confusingly: for PROcap data, we want the 5'
# ends of reads, but the reads appear reverse-complemented, so
# use the 3' end (and swap the files)
# (also for the G4RP paired-end reads)
if args.type in ['PROcap', 'paired']:
    flags = ['-3']
    # if this is paired-end reads, use the complete fragment
    if args.type in ['paired']:
        flags = ['-pc']
    # write out the individual strands (note that these are
    # swapped, compared to RNAseq and PROseq; see comment above)
    writeBigWig(args.BAM_input, args.bigWig_output + '.plus.bw',
        flags, normalize = not args.unnormalized,
        strand = '-', negate = False)
    # again, omitting negating
    writeBigWig(args.BAM_input, args.bigWig_output + '.minus.bw',
        flags, normalize = not args.unnormalized,
        strand = '+', negate = True)
    sys.exit(0)

