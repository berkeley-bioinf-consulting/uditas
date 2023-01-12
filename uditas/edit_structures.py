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
@click.option('-p', '--pos', help='guide cut position in amplicon', type=int, default=1000)
@click.option('-w', '--window', default=15, help='Window size for detecting indels', type=int)
@click.option('-n', '--num-output', help='Number of output structures (default: all structures)', type=int, default=None)
@click.option('-s', '--skip-umi', type=bool, is_flag=True, help='Skip UMI collapse')
@click.option('-u', '--umi-cluster', type=bool, is_flag=True, help='Print sequences clustered by UMI')
def structures(bam, umi, pos, window, num_output, skip_umi, umi_cluster):
    logger.info('Reading UMI barcodes')
    umi_dict = create_barcode_dict(umi)
    
    logger.info('Reading BAM file and creating pileups')
    samfile = pysam.AlignmentFile(bam, "rb")
    pu_start = samfile.pileup(contig='wt', start=pos-window, stop=pos+window, max_depth=100000, truncate=True, stepper='all')

    read_substr_dict = defaultdict(str)

    col = 0
    for pc in pu_start:
        col_seqs = {}
        max_len = 0
        for read in pc.pileups:
            qpos = read.query_position
            query = read.alignment.query_name
            umi_seq = umi_dict[query]
            if 'N' in umi_seq:
                # skip all reads with bad UMI
                continue
            query_id = (query, read.alignment.flag)
            #if query == 'M03561:432:GW22100230:1:2112:13907:13507':
            #    debug = True
            #else:
            debug = False
            if read.is_tail:
                # mark for deletion since it doesn't cover the full extent
                seq = '*'

            elif read.is_del:
                seq = '-'
            else:
                seq = read.alignment.query_sequence[qpos:qpos + read.indel + 1]
                max_len = max(max_len, len(seq))
            col_seqs[query_id] = seq
            if debug:
                sys.stderr.write(f'col: {col}, pos: {qpos}, del: {read.is_del}, indel: {read.indel}, col seq: "{seq}"\n')

        # first time through, create all the entries
        if col == 0:
            for query_id, col_seq in col_seqs.items():
                if col_seq == '*':
                    continue
                seq = read_substr_dict[query_id]
                seq += col_seq + ('.' * (max_len - len(col_seq)))
                read_substr_dict[query_id] = seq

        else:
            # subsequent iterations, only use existing
            # (prevent reads that start part-way through)
            to_del = []
            for query_id, seq in read_substr_dict.items():
                try:
                    col_seq = col_seqs[query_id]
                except KeyError:
                    # catch phantom deletions that don't have flag:
                    col_seq = '-'
                if col_seq == '*':
                    to_del.append(query_id)
                    continue
                # not sure how this continues to happen, but treat as del:
                if col_seq == '':
                    col_seq = '-'
                seq += col_seq + ('.' * (max_len - len(col_seq)))
                read_substr_dict[query_id] = seq
            for query_id in to_del:
                del(read_substr_dict[query_id])

        col += 1

    logger.info('Collapsing by UMI')
    umi_sets = defaultdict(Counter)
    for (query, _flag), seq in read_substr_dict.items():
        q_umi = umi_dict[query][0]
        umi_sets[q_umi][seq] += 1

    c = Counter()
    #umi_set = set()
    total_collapsed = 0
    for q_umi, umi_counter in umi_sets.items():
        for (seq, ct) in umi_counter.most_common(1):
#        if skip_umi:
#            total_collapsed += 1
#            c[seq] += 1
#        elif q_umi not in umi_set:
            #umi_set.add(q_umi)
            total_collapsed += 1
            c[seq] += 1

    print('Structure\tCount\tPct')
    all_seqs, all_counts = zip(*c.most_common(num_output))
    trimmed_seqs = trim_alignment(all_seqs)

    try:
        for (seq, ct) in zip(trimmed_seqs, all_counts):
            # all_seqs.append(seq)
            print(f'{seq}\t{ct}\t{ct/total_collapsed*100:.2f}')
    except:
        pass

    if umi_cluster:
        ct = 0
        for umi_seq, umi_counter in umi_sets.items():
            print('UMI ' + umi_seq)
            print('-' * 14)

            umi_seqs, umi_counts = zip(*umi_counter.most_common())
            trimmed_umi_seqs = trim_alignment(umi_seqs)
            umi_len = sum(umi_counts)
            for (seq, umi_ct) in zip(trimmed_umi_seqs, umi_counts):
                # all_seqs.append(seq)
                print(f'{seq}\t{umi_ct}\t{umi_ct/umi_len*100:.2f}')
            print()
            ct +=1
            #if ct == 20:
            #    break
    # print(trim_alignment(all_seqs[1:10]))
if __name__ == '__main__':
    logging.basicConfig(level=logging.INFO)
    structures()
