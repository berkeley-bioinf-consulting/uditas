import click
import random
import logging
import gzip
import os
import shutil
from typing import TextIO

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
@click.argument('source-dir', type=str)
@click.argument('target-dir', type=str)
@click.option('--make-cutadapt', help='Create cut adapt files', type=bool, default=False)
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


if __name__ == '__main__':
    cli()

