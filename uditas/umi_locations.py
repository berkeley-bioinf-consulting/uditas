import pysam
from collections import defaultdict, Counter
import gzip
import os
import click
import sys
import logging

logger = logging.getLogger()

def open_fastq_or_gz(filename):
    if filename.endswith(".fastq") and os.access(filename, os.F_OK):
        return open(filename, "r")
    elif filename.endswith(".fastq.gz") and os.access(filename, os.F_OK):
        return gzip.open(filename, "rt")
    elif filename.endswith(".fastq") and os.access(filename + ".gz", os.F_OK):
        return gzip.open(filename + ".gz", "rt")
    elif filename.endswith(".fastq.gz") and os.access(filename[:-3], os.F_OK):
        return open(filename[:-3], "r")
    raise IOError("Unknown file: " + filename)


def create_barcode_dict(filename): 
    barcode_file = open_fastq_or_gz(filename) 
    barcode_dict = dict()
    barcode_reads = zip(barcode_file)

    for header_barcode in barcode_reads:
        seq_barcode = next(barcode_reads)
        next(barcode_reads)
        qual_barcode = next(barcode_reads)
        barcode_dict[header_barcode[0].split()[0][1:]] = [seq_barcode[0].rstrip(), qual_barcode[0].rstrip()]

    return barcode_dict


def trim_alignment(seqs):
    seq_len = len(seqs[0])
    drop_pos = []

    for pos in range(seq_len):
        for seq in seqs:
            if pos == 0:
                assert len(seq) == seq_len, 'Sequences must be the same length'
            keep_pos = False
            if seq[pos] != '.':
                keep_pos = True
                break
        if not keep_pos:
            drop_pos.append(pos)
    last_pos = 0
    final_seqs = [''] * len(seqs)
    for pos in drop_pos:
        for i, seq in enumerate(seqs):
            final_seqs[i] += seq[last_pos:pos]
        last_pos = pos+1
    for i, seq in enumerate(seqs):
        final_seqs[i] += seq[last_pos:seq_len]
    return final_seqs


@click.command()
@click.argument('bam', type=str) # , help='BAM file for analysis')
@click.argument('umi', type=str) #, help='Sample UMI file')
def locations(bam, umi):
    logger.info('Reading UMI barcodes')
    umi_dict = create_barcode_dict(umi)
    
    logger.info('Reading BAM file and iterating')
    samfile = pysam.AlignmentFile(bam, "rb")
    seen_reads = set()
    umi_to_reads = defaultdict(Counter)

    for read in samfile.fetch():
        read_id = read.query_name 
        if read_id in seen_reads:
            continue
        seen_reads.add(read_id)
        read_umi = umi_dict[read_id][0]
        #logger.info(read_umi)
        umi_to_reads[read_umi][f'{read.reference_name}:{read.reference_start}'] += 1

    ct = 0
    for umi_seq, umi_counter in umi_to_reads.items():
        print('UMI ' + umi_seq)
        print('-' * 14)

        umi_seqs, umi_counts = zip(*umi_counter.most_common())
        #trimmed_umi_seqs = trim_alignment(umi_seqs)
        umi_len = sum(umi_counts)
        for (seq, umi_ct) in zip(umi_seqs, umi_counts):
            # all_seqs.append(seq)
            print(f'{seq}\t{umi_ct}\t{umi_ct/umi_len*100:.2f}')
        print()
        ct +=1


if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    locations()
