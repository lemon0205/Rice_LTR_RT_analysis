import os

pav_dir = 'PAV'
te_dir = 'TE-GFF'
output_dir = 'output'

bed_files = [f for f in os.listdir(pav_dir) if f.endswith('.final.bed')]
varieties = [f.replace('MSU_', '').replace('.final.bed', '') for f in bed_files]
outgroup = 'CG14'

# 长度相似阈值：允许的长度差异比例
LENGTH_SIMILARITY_RATIO = 0.3

def read_te_bed(file_path, family, source):
    tes_by_chrom = {}
    if not os.path.exists(file_path):
        return tes_by_chrom
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            te_id = parts[3]
            length = end - start
            if chrom not in tes_by_chrom:
                tes_by_chrom[chrom] = []
            tes_by_chrom[chrom].append({
                'start': start, 'end': end,
                'te_id': te_id, 'family': family, 'source': source,
                'length': length
            })
    for chrom in tes_by_chrom:
        tes_by_chrom[chrom].sort(key=lambda x: x['start'])
    return tes_by_chrom

def parse_variety_coords(coord_str):
    try:
        parts = coord_str.split(':')
        variety_chrom = parts[0]
        positions = parts[1]
        start_end = positions.split('-')
        start = int(start_end[0])
        end = int(start_end[1])
        chr_idx = variety_chrom.find('_Chr')
        if chr_idx == -1:
            return None, None, None, None
        variety = variety_chrom[:chr_idx]
        chrom = variety_chrom[chr_idx + 1:]
        return variety, chrom, start, end
    except:
        return None, None, None, None

# ========== Step 0: 读取所有数据 ==========
print("Step 0: 读取所有数据...")

# MSU参考基因组的TE
print("  读取MSU TE注释...")
msu_tes_list = []
for family in ['Copia', 'Gypsy']:
    bed_path = os.path.join(te_dir, family, f'{family}_bed', f'MSU_TE_{family}.bed')
    if not os.path.exists(bed_path):
        continue
    with open(bed_path, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            msu_tes_list.append({
                'te_id': parts[3],
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'family': family,
                'length': int(parts[2]) - int(parts[1])
            })
msu_tes_list.sort(key=lambda x: (x['chrom'], x['start']))
print(f"    MSU TE总数: {len(msu_tes_list)}")

# 各品种的TE注释
print("  读取各品种TE注释...")
variety_tes = {}
for variety in varieties:
    variety_tes[variety] = {}
    for family in ['Copia', 'Gypsy']:
        bed_path = os.path.join(te_dir, family, f'{family}_bed', f'{variety}_TE_{family}.bed')
        tes = read_te_bed(bed_path, family, variety)
        for chrom, te_list in tes.items():
            if chrom not in variety_tes[variety]:
                variety_tes[variety][chrom] = []
            variety_tes[variety][chrom].extend(te_list)
    total = sum(len(v) for v in variety_tes[variety].values())
    print(f"    {variety}: {total} TEs")

# 各品种的共线性文件 (coords.csv)
print("  读取各品种共线性文件...")
variety_coords = {}
for variety in varieties:
    coords_file = os.path.join(pav_dir, f'MSU_{variety}.bed.coords.csv')
    
    # 按MSU染色体分组存储共线性区间
    coord_intervals = []
    if os.path.exists(coords_file):
        with open(coords_file, 'r') as f:
            header = f.readline()
            for line in f:
                parts = line.strip().split(',')
                if len(parts) < 9:
                    continue
                ref_start = int(parts[0])
                ref_end = int(parts[1])
                query_start = int(parts[2])
                query_end = int(parts[3])
                ref_chrom = parts[6]
                query_contig = parts[7]
                
                coord_intervals.append({
                    'ref_chrom': ref_chrom,
                    'ref_start': ref_start,
                    'ref_end': ref_end,
                    'query_contig': query_contig,
                    'query_start': min(query_start, query_end),
                    'query_end': max(query_start, query_end),
                    'tag': parts[8],
                })
    
    coord_intervals.sort(key=lambda x: (x['ref_chrom'], x['ref_start']))
    variety_coords[variety] = coord_intervals
    print(f"    {variety}: {len(coord_intervals)} 个共线性区间")

# 各品种的PAV（用于辅助验证/覆盖）
print("  读取各品种PAV...")
variety_pav_index = {}
for variety in varieties:
    pav_file = os.path.join(pav_dir, f'MSU_{variety}.final.bed')
    intervals = []
    with open(pav_file, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            intervals.append({
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'pav_type': parts[5],
            })
    intervals.sort(key=lambda x: (x['chrom'], x['start']))
    variety_pav_index[variety] = intervals

# ========== Step 1: 共线性映射函数 ==========
def find_syntenic_region(coord_intervals, msu_chrom, msu_start, msu_end):
    """
    查找MSU坐标在品种中的共线性对应区域
    返回: list of (query_contig, query_start, query_end) 或空列表
    """
    regions = []
    for ci in coord_intervals:
        if ci['ref_chrom'] != msu_chrom:
            continue
        if ci['ref_start'] > msu_end:
            break
        
        # 检查是否有重叠
        if ci['ref_end'] >= msu_start and ci['ref_start'] <= msu_end:
            regions.append((ci['query_contig'], ci['query_start'], ci['query_end']))
    
    return regions

def find_matching_te_in_region(tes_by_chrom, regions, target_family, target_length, similarity_ratio=0.3):
    """
    在共线性区域内查找匹配的TE（同家族 + 相似长度）
    返回: 匹配到的TE列表 或 None(无匹配)
    """
    matched_tes = []
    min_len = target_length * (1 - similarity_ratio)
    max_len = target_length * (1 + similarity_ratio)
    
    for contig, reg_start, reg_end in regions:
        if contig not in tes_by_chrom:
            continue
        for te in tes_by_chrom[contig]:
            if te['family'] != target_family:
                continue
            if te['start'] > reg_end:
                break
            if te['end'] >= reg_start and te['start'] <= reg_end:
                if min_len <= te['length'] <= max_len:
                    matched_tes.append(te)
    
    return matched_tes if matched_tes else None

def query_pav_status(pav_intervals, te_chrom, te_start, te_end):
    """查询PAV状态"""
    for p in pav_intervals:
        if p['chrom'] != te_chrom:
            continue
        if p['start'] > te_end:
            break
        if p['end'] >= te_start and p['start'] <= te_end:
            return p['pav_type']
    return None

# ========== Step 2: 生成新的基因型矩阵 ==========
print("\n" + "="*60)
print("Step 2: 基于共线性映射生成TE基因型矩阵")
print("="*60)
print(f"  判断逻辑:")
print(f"    1. 先用PAV判断: INS→1, DEL→0")
print(f"    2. PAV未覆盖时用共线性映射:")
print(f"       找到共线性区域 → 在区域内查找同家族+相似长度TE")
print(f"       找到匹配TE → 1, 未找到 → 0")
print(f"       无共线性区域 → NA")

matrix_file = os.path.join(output_dir, 'te_genotype_matrix_final.csv')
with open(matrix_file, 'w') as f:
    headers = ['TE_ID', 'Family', 'Chromosome', 'Start', 'End', 'Length'] + varieties
    f.write(','.join(headers) + '\n')
    
    stats = {'pav_ins': 0, 'pav_del': 0, 'syn_match': 0, 'syn_nomatch': 0, 'na': 0}
    col_stats = {v: {'1': 0, '0': 0, 'NA': 0} for v in varieties}
    
    for i, te in enumerate(msu_tes_list):
        row = [te['te_id'], te['family'], te['chrom'], te['start'], te['end'], te['length']]
        
        for variety in varieties:
            # --- 第一步：先查PAV ---
            pav_type = query_pav_status(variety_pav_index[variety], te['chrom'], te['start'], te['end'])
            
            if pav_type == 'Insertion':
                row.append('1')
                col_stats[variety]['1'] += 1
                stats['pav_ins'] += 1
            elif pav_type == 'Deletion':
                row.append('0')
                col_stats[variety]['0'] += 1
                stats['pav_del'] += 1
            else:
                # --- 第二步：PAV未覆盖，用共线性映射 ---
                syn_regions = find_syntenic_region(
                    variety_coords[variety], te['chrom'], te['start'], te['end']
                )
                
                if not syn_regions:
                    # 无共线性区域 → NA
                    row.append('NA')
                    col_stats[variety]['NA'] += 1
                    stats['na'] += 1
                else:
                    # 有共线性区域，查找匹配的TE
                    matched = find_matching_te_in_region(
                        variety_tes.get(variety, {}),
                        syn_regions,
                        te['family'],
                        te['length'],
                        LENGTH_SIMILARITY_RATIO
                    )
                    
                    if matched:
                        row.append('1')
                        col_stats[variety]['1'] += 1
                        stats['syn_match'] += 1
                    else:
                        row.append('0')
                        col_stats[variety]['0'] += 1
                        stats['syn_nomatch'] += 1
        
        f.write(','.join(str(x) for x in row) + '\n')
        
        if (i + 1) % 10000 == 0:
            print(f"    已处理 {i+1}/{len(msu_tes_list)} 个TE")

print(f"\n  保存到 {matrix_file}")
print(f"  总行数(TE数): {len(msu_tes_list)}")

print(f"\n  全局统计:")
print(f"    PAV-Insertion→1: {stats['pav_ins']}")
print(f"    PAV-Deletion→0:  {stats['pav_del']}")
print(f"    共线性-匹配→1:   {stats['syn_match']}")
print(f"    共线性-未匹配→0: {stats['syn_nomatch']}")
print(f"    无共线性→NA:     {stats['na']}")

print(f"\n  各品种列统计:")
print(f"    {'品种':<12} {'存在(1)':>10} {'缺失(0)':>10} {'NA':>8}")
for v in varieties:
    cs = col_stats[v]
    total = cs['1'] + cs['0'] + cs['NA']
    print(f"    {v:<12} {cs['1']:>10} {cs['0']:>10} {cs['NA']:>8}")

# ========== Step 3: dTE分析 ==========
print("\n" + "="*60)
print("Step 3: dTE分析（基于外群CG14）")
print("="*60)

dte_file = os.path.join(output_dir, 'te_dte_analysis_final.csv')

with open(matrix_file, 'r') as fin:
    header = fin.readline().strip().split(',')
    
    col_idx = {v: header.index(v) for v in varieties}
    cg14_col = col_idx[outgroup]
    
    dte_stats = {v: {'Gain': 0, 'Loss': 0, 'Ancestral': 0, 'NA': 0} for v in varieties if v != outgroup}
    
    with open(dte_file, 'w') as fdte:
        dte_headers = ['TE_ID', 'Family', 'Chromosome', 'Start', 'End', outgroup] + \
                      [v for v in varieties if v != outgroup]
        fdte.write(','.join(dte_headers) + '\n')
        
        for line in fin:
            parts = line.strip().split(',')
            if len(parts) < len(header):
                continue
            
            te_info = parts[:6]
            
            cg14_val_str = parts[cg14_col]
            
            row_dte = list(te_info) + [cg14_val_str]
            
            for variety in varieties:
                if variety == outgroup:
                    continue
                
                var_val_str = parts[col_idx[variety]]
                
                if cg14_val_str == 'NA' or var_val_str == 'NA':
                    dte = 'NA'
                    dte_stats[variety]['NA'] += 1
                elif cg14_val_str == '0' and var_val_str == '1':
                    dte = 'Gain'
                    dte_stats[variety]['Gain'] += 1
                elif cg14_val_str == '1' and var_val_str == '0':
                    dte = 'Loss'
                    dte_stats[variety]['Loss'] += 1
                else:
                    dte = 'Ancestral'
                    dte_stats[variety]['Ancestral'] += 1
                
                row_dte.append(dte)
            
            fdte.write(','.join(row_dte) + '\n')

print(f"  保存到 {dte_file}")

# ========== Step 4: 统计汇总 ==========
print("\n" + "="*60)
print("Step 4: 统计汇总")
print("="*60)

stats_file = os.path.join(output_dir, 'te_dte_stats_final.csv')
with open(stats_file, 'w') as f:
    f.write('Variety,Gain,Loss,Ancestral,NA,Total\n')
    for variety in sorted(dte_stats.keys()):
        s = dte_stats[variety]
        total = s['Gain'] + s['Loss'] + s['Ancestral'] + s['NA']
        f.write(f"{variety},{s['Gain']},{s['Loss']},{s['Ancestral']},{s['NA']},{total}\n")

print(f"\n{'Variety':<12} {'Gain':>8} {'Loss':>8} {'Ancestral':>10} {'NA':>6} {'Total':>8}")
for variety in sorted(dte_stats.keys()):
    s = dte_stats[variety]
    total = s['Gain'] + s['Loss'] + s['Ancestral'] + s['NA']
    print(f"{variety:<12} {s['Gain']:>8} {s['Loss']:>8} {s['Ancestral']:>10} {s['NA']:>6} {total:>8}")

print(f"\n{'='*60}")
print("分析完成！生成的文件:")
print(f"{'='*60}")
print(f"  output/")
print(f"  ├── te_genotype_matrix_final.csv  ← 新版矩阵（PAV + 共线性映射）")
print(f"  ├── te_dte_analysis_final.csv     ← dTE分析结果")
print(f"  └── te_dte_stats_final.csv        ← dTE统计结果")
