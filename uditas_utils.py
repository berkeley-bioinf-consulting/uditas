import click
import random
import logging
import gzip
import os
import shutil
from typing import TextIO
from collections import defaultdict


import pandas as pd

BASES = 'ACGT'

logger = logging.getLogger(__name__)


@click.group()
@click.option('--debug/--no-debug', default=False)
def cli(debug):
    logging.basicConfig(level=logging.DEBUG if debug else logging.INFO)


def create_umi_file(fq: str, output: str, length: int, qual_value: str):
    with open_gz_or_txt(fq, 'r') as fq_file, open_gz_or_txt(output, 'w') as umi_file:
        write_umi(fq_file, umi_file, length, qual_value)


def write_umi(fq_file: TextIO, umi_file: TextIO, length: int, qual_value: str):
    umi_set = set()
    for header in fq_file:
        umi, qual = make_umi(length, qual_value)
        while umi in umi_set:
            umi, qual = make_umi(length, qual_value)
        umi_set.add(umi)
        umi_file.write('\n'.join([header.strip(), umi, '+', qual, '']))
        next(fq_file)
        next(fq_file)
        next(fq_file)


def make_umi(length, qual_value):
    return ''.join([random.choice(BASES) for x in range(length)]), qual_value * length


def open_gz_or_txt(f, mode):
    if f.endswith('.gz'):
        return gzip.open(f, f'{mode}t')
    else:
        return open(f, mode)


def read_samples(infile):
    samples = pd.read_csv(infile)
    samples['R1'] = samples['Sample'] + '_R1_001.fastq.gz'
    samples['R2'] = samples['Sample'] + '_R2_001.fastq.gz'
    samples['sample_dir'] = samples['index_I1'] + '_' + samples['index_I2']
    return samples


def make_sample_dirs(samples, base_dir, make_cut):
    for sample_dir in samples['sample_dir']:
        os.makedirs(os.path.join(base_dir, sample_dir), exist_ok=True)
        os.makedirs(os.path.join(base_dir, sample_dir, 'fastq_files'), exist_ok=True)
        if make_cut:
            os.makedirs(os.path.join(base_dir, sample_dir, 'cutadapt_files'), exist_ok=True)


def move_sample_files(samples, source_dir, base_dir, make_cut):
    for index, row in samples.iterrows():
        r1_dst = os.path.join(base_dir, row['sample_dir'], 'fastq_files', row['sample_dir'] + '_R1.fastq.gz')
        r2_dst = os.path.join(base_dir, row['sample_dir'], 'fastq_files', row['sample_dir'] + '_R2.fastq.gz')
        shutil.copyfile(os.path.join(source_dir, row['R1']), r1_dst)
        shutil.copyfile(os.path.join(source_dir, row['R2']), r2_dst)
        if make_cut:
            shutil.copyfile(r1_dst, os.path.join(base_dir, row['sample_dir'], 'cutadapt_files', row['sample_dir'] + '_R1.trimmed.fastq.gz'))
            shutil.copyfile(r2_dst, os.path.join(base_dir, row['sample_dir'], 'cutadapt_files', row['sample_dir'] + '_R2.trimmed.fastq.gz'))


def make_umi_files(samples, base_dir):
    for s in samples['sample_dir']:
        create_umi_file(os.path.join(base_dir, s, 'fastq_files', s + '_R1.fastq.gz'), os.path.join(base_dir, s, 'fastq_files', s + '_umi.fastq.gz'), 16, 'F')


@click.command('preprocess')
@click.argument('sample-info', type=str)
@click.argument('source-dir', type=str)
@click.argument('target-dir', type=str)
@click.option('--make-cutadapt', help='Assume fq files are trimmed and copy to cut-adapt directory', type=bool, default=False)
def preprocess(sample_info: str, source_dir: str, target_dir: str, make_cutadapt: bool):
    if not os.path.exists(source_dir):
        raise ValueError('Source dir must exist')
    sample_info = os.path.join(target_dir, 'sample_info.csv')
    if not os.path.exists(sample_info):
        raise FileNotFoundError(f'File {sample_info} must be present')

    samples = read_samples(sample_info)
    make_sample_dirs(samples, target_dir, make_cutadapt)
    move_sample_files(samples, source_dir, target_dir, make_cutadapt)
    make_umi_files(samples, target_dir)


cli.add_command(preprocess)

@click.command('demux')
@click.argument('sample-info', type=str)
@click.argument('source-dir', type=str)
@click.argument('target-dir', type=str)
@click.option('--make-cutadapt', help='Assume fq files are trimmed and copy to cut-adapt directory', type=bool, default=False)
def demux(sample_info: str, source_dir: str, target_dir: str, make_cutadapt: bool):
    if not os.path.exists(source_dir):
        raise ValueError('Source dir must exist')
    sample_info = os.path.join(target_dir, 'sample_info.csv')
    if not os.path.exists(sample_info):
        raise FileNotFoundError(f'File {sample_info} must be present')

    samples = read_samples(sample_info)
    make_sample_dirs(samples, target_dir, make_cutadapt)

    
cli.add_command(demux)

if __name__ == '__main__':
    cli()

file_i1 = 'Un_L001_I1_001.fastq.gz'
file_i2 = 'Un_L001_I2_001.fastq.gz'
file_r1 = 'Un_L001_R1_001.fastq.gz'
file_r2 = 'Un_L001_R2_001.fastq.gz'
file_r3 = 'Un_L001_R3_001.fastq.gz'

def create_demux_handles(samples, infile_pat='Un_L001_{read_type}_001.fastq.gz', outfile_pat='{sample}_L001_{read_type}_001.fastq.gz'):
    read_types = ['I1', 'I2', 'R1', 'R2', 'R3']


    handles = defaultdict(dict)
    for sample in samples['sample']:
        # sample = dct['sample']
        for read_type in read_types:
            # put these into the correct directories.
            handles[sample][read_type] = gzip.open(outfile_pat.format(sample=sample, read_type=read_type), 'wt')  # f'{sample}_L001_{read_type}_001.fastq.gz', 'wt')

    for read_type in read_types:
        handles['Un'][read_type] = gzip.open(outfile_pat.format(sample='Un', read_type=read_type), 'wt')  # f'Un2_L001_{read_type}_001.fastq.gz', 'wt')

    in_handles = {rt: gzip.open(infile_pat.format(read_type=rt,'rt') for rt in read_types}
    return in_handles, handles

files = {
        'I1': file_i1,
        'I2': file_i2,
        'R1': file_r1,
        'R2': file_r2,
        'R3': file_r3
}

sample_file = 'samples.csv'
samples = [
        {'index':'CTCTCTAT','sample':'S1'},
        {'index':'CTAAGCCT','sample':'S2'},
        {'index':'ACTGCATA','sample':'S3'},
        {'index':'ACTGATGG','sample':'S4'},
        {'index':'AAGGAGTA','sample':'S5'},
        {'index':'GTACCTAG','sample':'S6'},
        {'index':'TATCCTCT','sample':'S7'},
        {'index':'AGAGTAGA','sample':'S9'},
        {'index':'GTAAGGAG','sample':'S9_2'},
        {'index':'GACATTGT','sample':'S8_maybe'}
]


samp_dict = {}

index = dct['index']

samp_dict[index] = sample

line_ct = 0
total_ct = 0
buf = defaultdict(str)
while True:
    line_ct += 1
    for rt, fh in in_handles.items():
        line = fh.readline()
        if not line:
            break

        buf[rt] += line

        if rt == 'I2' and line_ct == 2:
            index = line.strip()
            try:
                sample = samp_dict[index]
            except KeyError:
                sample = 'Un'
    if line_ct == 4:
        line_ct = 0

        fh_list = handles[sample]
        for rt, fh in fh_list.items():
            fh.write(buf[rt])
        buf = defaultdict(str)

    total_ct += 1
    if total_ct % 100000 == 0:
        print(f'read {total_ct} mates')
    if not line:
        break

for fh in in_handles.values():
    fh.close()

for dct in handles:
    for fh in dct.values():
        fh.close()

dna_tr = str.maketrans('ACGT', 'TGCA')

def rc(seq):
    return seq.translate(dna_tr)[::-1]

sample_file = 'samples.csv'
samples = [
        {'index': 'GCCTTA', 'index2':'CTCTCTAT','sample':'S1', 'out_name': 'A1-I701_A01'},
        {'index': 'CTGCCT', 'index2':'CTAAGCCT','sample':'S2', 'out_name': 'A2-I701_A02'},
        {'index': 'GAGTCC', 'index2':'ACTGCATA','sample':'S3', 'out_name': 'A3-I701_A03'}, # sample 3
        {'index': 'AGAGAG', 'index2':'ACTGATGG','sample':'S4', 'out_name': 'A4-I701_A04'},
        {'index': 'TCAGGA', 'index2':'AAGGAGTA','sample':'S5', 'out_name': 'A5-I701_A05'},
        {'index': 'TCTCTG', 'index2':'GTACCTAG','sample':'S6', 'out_name': 'A6-I701_A06'},
        {'index': 'AGTACG', 'index2':'TATCCTCT','sample':'S7', 'out_name': 'A7-I701_A07'},
        {'index': 'TGCCTA', 'index2':'GACATTGT','sample':'S8_maybe', 'out_name': 'A8-I701_A08'},  # sample 8?
        {'index': 'AGTACG', 'index2':'AGAGTAGA','sample':'S9', 'out_name': 'A9-I701_A09'},
        {'index': 'AGTACG', 'index2':'GTAAGGAG','sample':'S9_2', 'out_name': 'A10-I701_A10'}
]

read_types = ['I1', 'I2', 'R1', 'R2', 'R3']

samp_dict = {}
samp2_dict = {}

handles = defaultdict(lambda: defaultdict(dict))
for dct in samples:
    index = rc(dct['index'])
    index2 = dct['index2']
    sample = dct['sample']
    out = dct['out_name']
    samp_dict[index] = sample
    print(f'found sample {sample} with indexes {index} and {index2}.')
    samp2_dict[index2] = sample
    try:
        os.mkdir(out)
    except:
        pass

    for read_type in read_types:
        handles[sample]['in'][read_type] = gzip.open(f'{sample}_L001_{read_type}_001.fastq.gz', 'rt')
        out_type = read_type
        if out_type == 'R2':
            out_type = 'umi'
        elif out_type == 'R3':
            out_type = 'R2'
        handles[sample]['out'][read_type] = gzip.open(f'{out}/{out}_{out_type}.fastq.gz', 'wt')
        handles[sample]['Un'][read_type] = gzip.open(f'{out}/Un_{out}_{read_type}.fastq.gz', 'wt')

for sample in handles:
    line_ct = 0
    total_ct = 0
    buf = defaultdict(str)
    in_handles = handles[sample]['in']
    out_handles = handles[sample]['out']
    un_handles = handles[sample]['Un']
    samp_i1 = 'Un'
    samp_i2 = 'Un'
    while True:
        line_ct += 1
        for rt, fh in in_handles.items():
            line = fh.readline()
            if not line:
                break

            buf[rt] += line

            if line_ct == 2:
                index = line.strip()
                if rt == 'I1':
                    try:
                        samp_i1 = samp_dict[index]
                    except KeyError:
                        samp_i1 = 'Un'
                elif rt == 'I2':
                    try:
                        samp_i2 = samp2_dict[index]
                    except KeyError:
                        samp_i2 = 'Un'
        if line_ct == 4:
            line_ct = 0
            if samp_i1 == sample and samp_i2 == sample:
                fh_list = out_handles
            else:
                fh_list = un_handles
            for rt, fh in fh_list.items():
                fh.write(buf[rt])
            buf = defaultdict(str)
            samp_i1 = 'Un'
            samp_i2 = 'Un'

        total_ct += 1
        if total_ct % 100000 == 0:
            print(f'read {total_ct} mates')
        if not line:
            break

    for rt, fh in in_handles.items():
        fh.close()
    for rt, fh in out_handles.items():
        fh.close()
    for rt, fh in un_handles.items():
        fh.close()
    print(f'Sample {sample} complete')
