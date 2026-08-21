#!/usr/bin/env python3
"""
HDF5 ファイル群のカバレッジ一覧を出力する。
出力は mean_depth 昇順のタブ区切りテキスト (stdout)。

使い方:
  python hdf5_depths.py /path/to/hdf5_dir [--sort filename|samplename|depth|zero]
"""

import argparse
import sys
import numpy as np
import h5py
from pathlib import Path


def read_stats(path):
    with h5py.File(path, 'r') as f:
        raw = f['/sample_metadata/sample_name'][()]
        name = raw[0].decode() if hasattr(raw[0], 'decode') else str(raw[0])

        if '/counts/values' in f:
            counts = f['/counts/values'][()].astype(float).ravel()
        elif '/counts/counts' in f:
            counts = f['/counts/counts'][()].astype(float).ravel()
        else:
            raise KeyError('counts データセットが見つかりません')

    return name, counts.mean(), (counts == 0).mean()


def main():
    parser = argparse.ArgumentParser(description='HDF5 depth 一覧出力')
    parser.add_argument('hdf5_dir', help='HDF5 ファイルが入っているディレクトリ')
    parser.add_argument('--sort', choices=['filename', 'samplename', 'depth', 'zero'],
                        default='depth', help='ソートキー (デフォルト: depth)')
    args = parser.parse_args()

    files = sorted(Path(args.hdf5_dir).glob('**/*.hdf5'))
    if not files:
        print('ERROR: .hdf5 ファイルが見つかりません', file=sys.stderr)
        sys.exit(1)

    rows = []
    for f in files:
        try:
            name, depth, zero_frac = read_stats(f)
            rows.append((f.name, name, depth, zero_frac))
        except Exception as e:
            print(f'SKIP {f.name}: {e}', file=sys.stderr)

    sort_key = {'filename': 0, 'samplename': 1, 'depth': 2, 'zero': 3}[args.sort]
    rows.sort(key=lambda r: r[sort_key])

    depths = [r[2] for r in rows]
    mean_d = np.mean(depths)
    std_d  = np.std(depths)
    med_d  = np.median(depths)

    #print(f'# サンプル数={len(rows)}  mean={mean_d:.1f}  median={med_d:.1f}  std={std_d:.1f}  sort={args.sort}')
    print('\t'.join(['filename', 'sample_name', 'mean_depth', 'zero_pct']))
    for fname, sname, depth, zfrac in rows:
        print(f'{fname}\t{sname}\t{depth:.2f}\t{zfrac*100:.1f}')


if __name__ == '__main__':
    main()
