# -*- coding: utf-8 -*-
"""
SELEX 文库序列分析完整流程（有 Tag / 无 Tag 自动适配）

模式判断规则：
  仅根据配置区 second_half_fields 是否填写 Tag：
    - second_half_fields 有内容 → 有 Tag 模式，按 Tag 分 library_1 ~ library_N
    - second_half_fields = []    → 无 Tag 模式，全部序列归入 library_1

注意：
  adapter、前引物、后引物、PCR 引物、random_length 等参数全部由用户手动填写。
  步骤2：count≥MIN_RETENTION_COUNT 过滤 → 各文库合并 → cd-hit-est 聚类 →
          输出：按家族丰度排序的全部序列 fasta、家族代表序列 fasta、家族统计 CSV
  步骤3：family_statistics.csv → dna_families_with_seq.txt → 渐进式多序列比对 → sequence logo (png+pdf)

使用方法：
  只需修改下方【配置区】中的测序数据路径、接头、前后引物、文库识别字段(tag)、
  软件路径等参数，然后直接运行本脚本即可跑完全流程。
"""

import os
import subprocess
import csv
import warnings
from collections import defaultdict, Counter

import regex as re     # 支持模糊匹配/错配正则: pip install regex
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import logomaker
import difflib

try:
    from Bio import SeqIO  # 需要安装 biopython: pip install biopython
except ImportError:
    SeqIO = None

warnings.filterwarnings('ignore')

# ================================【Web/API 配置读取】================================
# 本版本由网页/API 自动生成 config.json，并通过 SELEX_CONFIG_PATH 注入。
# 本地直接运行时也可以设置 SELEX_CONFIG_PATH 指向一个 JSON 配置文件。

import json
import shutil

CONFIG_PATH = os.environ.get("SELEX_CONFIG_PATH", "config.json")

def load_config():
    with open(CONFIG_PATH, "r", encoding="utf-8") as f:
        return json.load(f)

CFG = load_config()

# 1. 测序数据与输出路径
input_file1 = CFG["input_file1"]
input_file2 = CFG["input_file2"]
final_output_dir = CFG["final_output_dir"]
OUTPUT_DIR = CFG["OUTPUT_DIR"]
output_prefix = CFG["output_prefix"]

# 2. 接头 / 前后引物 / Tag
adapter1 = CFG["adapter1"]
adapter2 = CFG["adapter2"]
extract_begin_seq = CFG["extract_begin_seq"]
extract_end_seq = CFG["extract_end_seq"]
random_length = int(CFG["random_length"])
second_half_fields = CFG.get("second_half_fields", [])
reverse_primer = CFG.get("reverse_primer", "")
REMOVE_PCR_PRIMER = bool(CFG.get("REMOVE_PCR_PRIMER", True))
HAS_TAG = len(second_half_fields) > 0

# 3. 提取设置
max_mismatch = int(CFG.get("max_mismatch", 2))
FLASH_MAX_OVERLAP = str(CFG.get("FLASH_MAX_OVERLAP", "250"))

# 4. CD-HIT
SIM_THRESHOLD = float(CFG.get("SIM_THRESHOLD", 0.80))
THREADS = int(CFG.get("THREADS", 8))
MIN_RETENTION_COUNT = int(CFG.get("MIN_RETENTION_COUNT", 3))
select_file_names = CFG.get("select_file_names", None)

# 5. motif
MEMBER_COUNT_THRESH = int(CFG.get("MEMBER_COUNT_THRESH", 50))
SINGLE_SEQ_COUNT_THRESH = int(CFG.get("SINGLE_SEQ_COUNT_THRESH", 500))
MIN_COL_COVERAGE = float(CFG.get("MIN_COL_COVERAGE", 0.3))

# 6. 软件路径：优先使用配置，否则从 PATH 查找
CUTADAPT_PATH = CFG.get("CUTADAPT_PATH") or shutil.which("cutadapt") or "cutadapt"
FLASH_PATH = CFG.get("FLASH_PATH") or shutil.which("flash") or "flash"
CDHIT_PATH = CFG.get("CDHIT_PATH") or shutil.which("cd-hit-est") or "cd-hit-est"

# ================================【以下为流程代码，一般无需修改】===============================


DNA_COLORS = {'A': '#CC0000', 'T': '#00A000', 'C': '#0000CC', 'G': '#FFB300'}
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'Arial', 'Helvetica']

# ----------------- 通用工具函数 -----------------
def reverse_complement(sequence):
    """计算 DNA 序列的反向互补"""
    complement = str.maketrans('ATCGNatcgn', 'TAGCNtagcn')
    return sequence.translate(complement)[::-1]

def run_command(cmd, desc):
    """系统命令执行包装"""
    print(f"📌 {desc}...")
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True, errors='ignore')
        return True
    except Exception as e:
        print(f"❌ 命令执行失败: {desc} | 报错: {str(e)}")
        return False

def compile_primer_regex(primer, max_err):
    """构造带错配容忍的正则模式 (使用 regex 库)"""
    return re.compile(f"({primer}){{e<={max_err}}}")

def read_fasta_native(fa_path):
    """读取fasta，返回dict{seq_id:sequence}"""
    seq_dict = {}
    if not os.path.exists(fa_path):
        raise FileNotFoundError(f"文件不存在：{fa_path}")
    with open(fa_path, "r", encoding="utf-8") as f:
        lines = f.read().splitlines()
    sid = None
    buf = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if sid is not None and len(buf) > 0:
                seq_dict[sid] = "".join(buf)
            sid = line[1:]
            buf = []
        else:
            buf.append(line)
    if sid is not None and len(buf) > 0:
        seq_dict[sid] = "".join(buf)
    return seq_dict

def extract_count_from_id(seq_id):
    """从序列ID中提取count数字，兼容两种fasta头格式：
    >Seq_1_reads_828  和  >Seq_1_counts_828
    提取失败返回None"""
    for marker in ("_counts_", "_reads_"):
        if marker in seq_id:
            try:
                return int(seq_id.split(marker)[-1])
            except ValueError:
                return None
    return None


# ----------------- 步骤1a：文库字段归类 -----------------
def find_library_field(seq, primer_spans):
    """在整条序列中确定文库归属（返回 1~N，找不到返回 0）。

    不依赖文库字段必须位于序列两端（真实文库中字段可能被 PCR 引物
    或其他序列包在中间）。查找顺序：
      1) 左侧引物（begin/rc_begin）上游 40 bp 窗口：字段的 RC 形式优先，
         再正向形式；
      2) 右侧引物（end/rc_end）下游 40 bp 窗口：字段本身优先，再 RC 形式；
      3) 端点 30 bp 窗口兜底；
      4) 全序列兜底：仅排除引物匹配区间内部（不向外扩展），避免引物
         内部子串（如 end 引物 GATCTACCTTCGGGTTTGGC 含 TCTAC）误判文库。
    """
    # 无 Tag：不做 Tag 搜索，所有成功提取的序列统一归入 library_1
    if not HAS_TAG:
        return 1

    n_libs = len(second_half_fields)
    first_fields = [reverse_complement(f) for f in second_half_fields]
    L = len(seq)

    if primer_spans:
        sorted_spans = sorted(primer_spans, key=lambda x: x[0])
        left_start = sorted_spans[0][0]
        right_end = sorted_spans[-1][1]
        # 1) 左侧引物上游 40 bp：字段 RC 形式优先，再正向形式
        region = seq[max(0, left_start - 40):left_start]
        for i, f in enumerate(first_fields, 1):
            if f in region:
                return i
        for i, f in enumerate(second_half_fields, 1):
            if f in region:
                return i
        # 2) 右侧引物下游 40 bp：字段本身优先，再 RC 形式
        region = seq[right_end:right_end + 40]
        for i, f in enumerate(second_half_fields, 1):
            if f in region:
                return i
        for i, f in enumerate(first_fields, 1):
            if f in region:
                return i

    # 3) 端点 30 bp 窗口兜底
    for i, f in enumerate(first_fields, 1):
        if f in seq[:30]:
            return i
    for i, f in enumerate(second_half_fields, 1):
        if f in seq[-30:]:
            return i

    # 4) 全序列兜底（仅排除引物区间内部）
    excluded = set()
    for s, e in primer_spans:
        excluded.update(range(max(0, s), min(L, e)))
    for i, f in enumerate(second_half_fields, 1):
        pos = seq.find(f)
        if pos != -1 and not any((pos + j) in excluded for j in range(len(f))):
            return i
    for i, f in enumerate(first_fields, 1):
        pos = seq.find(f)
        if pos != -1 and not any((pos + j) in excluded for j in range(len(f))):
            return i
    return 0


# ----------------- 步骤1b：核心提取与按文库归类输出 -----------------
def extract_n45_and_classify(fastq_file, out_dir):
    """按用户设定的 random_length 提取随机区；有 Tag 按 Tag 分类，无 Tag 全部归入 library_1。"""
    os.makedirs(out_dir, exist_ok=True)

    # 有 Tag：N 个 Tag 对应 N 个 library
    # 无 Tag：统一使用 library_1
    n_libs = len(second_half_fields) if HAS_TAG else 1

    txt_handles = {}
    for i in range(1, n_libs + 1):
        txt_handles[i] = open(os.path.join(out_dir, f"library_{i}_cleaned_sequences.txt"), "w")
    counts = {i: defaultdict(int) for i in range(1, n_libs + 1)}

    # 构建正向与反向互补引物的模糊匹配器
    fwd_begin = compile_primer_regex(extract_begin_seq, max_mismatch)
    fwd_end   = compile_primer_regex(extract_end_seq, max_mismatch)
    rc_begin  = compile_primer_regex(reverse_complement(extract_end_seq), max_mismatch)
    rc_end    = compile_primer_regex(reverse_complement(extract_begin_seq), max_mismatch)

    total_records = 0
    extracted_count = 0
    unclassified_count = 0

    if HAS_TAG:
        print(f"🔄 开始逐条提取 N{random_length} 核心序列并按 {len(second_half_fields)} 个 Tag 归类...")
    else:
        print(f"🔄 开始逐条提取 N{random_length} 核心序列（无 Tag：全部归入 library_1）...")

    with open(fastq_file, "r") as in_f:
        while True:
            lines = [in_f.readline() for _ in range(4)]
            if not lines[3]:
                break

            total_records += 1
            seq = lines[1].strip().upper()
            extracted_seq = ""
            primer_spans = []

            # 方案 A: 正向匹配 (5' 前引 -> 3' 后引)
            m_begin = fwd_begin.search(seq)
            m_end   = fwd_end.search(seq)
            if m_begin:
                primer_spans.append(m_begin.span())
            if m_end:
                primer_spans.append(m_end.span())
            if m_begin and m_end and m_begin.end() < m_end.start():
                candidate = seq[m_begin.end():m_end.start()]
                if abs(len(candidate) - random_length) <= 3:
                    extracted_seq = candidate

            # 方案 B: 双向文库中的反向互补链匹配
            if not extracted_seq:
                m_rc_b = rc_begin.search(seq)
                m_rc_e = rc_end.search(seq)
                if m_rc_b:
                    primer_spans.append(m_rc_b.span())
                if m_rc_e:
                    primer_spans.append(m_rc_e.span())
                if m_rc_b and m_rc_e and m_rc_b.end() < m_rc_e.start():
                    candidate = seq[m_rc_b.end():m_rc_e.start()]
                    if abs(len(candidate) - random_length) <= 3:
                        extracted_seq = reverse_complement(candidate)

            # 方案 C (兜底策略): 仅匹配单侧前引并截取 45 bp
            if not extracted_seq and m_begin:
                start_p = m_begin.end()
                if start_p + random_length <= len(seq):
                    extracted_seq = seq[start_p : start_p + random_length]

            if not extracted_seq:
                continue
            extracted_count += 1

            # 按文库字段归类并写入对应文库
            lib = find_library_field(seq, primer_spans)
            if lib == 0:
                unclassified_count += 1
                continue
            txt_handles[lib].write(f"{extracted_seq}\n")
            counts[lib][extracted_seq] += 1

    for i in range(1, n_libs + 1):
        txt_handles[i].close()
        sorted_seqs = sorted(counts[i].items(), key=lambda x: x[1], reverse=True)
        with open(os.path.join(out_dir, f"library_{i}_dedup_sorted.fasta"), "w") as fasta_f:
            for rank, (s, c) in enumerate(sorted_seqs, 1):
                fasta_f.write(f">Seq_{rank}_reads_{c}\n{s}\n")
        print(f"✅ 文库 {i} 输出完成: 提取 {sum(counts[i].values())} 条 | unique {len(counts[i])} 条 | 已输出至 {out_dir}")

    extraction_rate = extracted_count / total_records if total_records else 0
    unclassified_rate = unclassified_count / extracted_count if extracted_count else 0
    print(f"✅ 提取完成！总输入 {total_records} | 成功提取 N{random_length} {extracted_count} "
          f"(有效率 {extraction_rate:.2%}) | 未归类 {unclassified_count} "
          f"({unclassified_rate:.2%})")


# ----------------- 步骤1c：自动化预处理主流程 -----------------
def process_pipeline():
    os.makedirs(final_output_dir, exist_ok=True)
    os.makedirs(os.path.dirname(output_prefix), exist_ok=True)

    # 1. QC & Cutadapt 去接头
    run_command([CUTADAPT_PATH, "-a", adapter1, "-q", "20", "-o", f"{output_prefix}_R1_trimmed.fq", input_file1], "R1 过滤与去接头")
    run_command([CUTADAPT_PATH, "-a", adapter2, "-q", "20", "-o", f"{output_prefix}_R2_trimmed.fq", input_file2], "R2 过滤与去接头")

    # 2. FLASH 合并双端 Reads
    run_command([
        FLASH_PATH,
        f"{output_prefix}_R1_trimmed.fq",
        f"{output_prefix}_R2_trimmed.fq",
        "-o", f"{output_prefix}_flash",
        "-M", FLASH_MAX_OVERLAP
    ], "FLASH 拼接双端序列")

    merged_fastq = f"{output_prefix}_flash.extendedFrags.fastq"
    if not os.path.exists(merged_fastq):
        print("❌ FLASH 合并文件不存在，请检查 FLASH 执行日志。")
        return

    # 3. 外侧 PCR 引物处理
    # 是否执行完全由 REMOVE_PCR_PRIMER 控制，与有无 Tag 无关。
    if REMOVE_PCR_PRIMER and reverse_primer.strip():
        primer1 = reverse_complement(reverse_primer)
        run_command([
            CUTADAPT_PATH, "-g", primer1,
            "-o", f"{output_prefix}_primers_trimmed.fq",
            merged_fastq
        ], "去除正向 PCR 引物")

        run_command([
            CUTADAPT_PATH, "-a", reverse_primer,
            "-o", f"{output_prefix}_final_pre_split.fq",
            f"{output_prefix}_primers_trimmed.fq"
        ], "去除反向 PCR 引物")

        split_input = f"{output_prefix}_final_pre_split.fq"
        if not os.path.exists(split_input):
            print("⚠️ 去 PCR 引物后文件不存在，改用 FLASH 合并文件直接提取。")
            split_input = merged_fastq
    else:
        split_input = merged_fastq
        if not reverse_primer.strip():
            print("ℹ️ 未填写 reverse_primer：跳过去除外侧 PCR 引物。")
        else:
            print("ℹ️ REMOVE_PCR_PRIMER=False：跳过去除外侧 PCR 引物。")

    # 4. 提取随机区；random_length 完全使用用户配置值
    extract_n45_and_classify(split_input, final_output_dir)


# ----------------- 步骤2a：预处理过滤 fasta -----------------
def filter_fasta_by_count(input_fa, output_fa, min_count=3):
    """
    预处理：过滤fasta，丢弃count <= min_count-1的序列
    fasta ID格式 >xxx_reads_5 或 >xxx_counts_5
    同时给每个序列ID前面加上来源前缀（从文件名提取，如 library_1）
    :param input_fa:原始输入fasta
    :param output_fa:过滤之后输出fasta
    :param min_count:最小保留count，保留count>=min_count
    """
    # 从文件名提取来源前缀，例如 library_1_dedup_sorted.fasta -> library_1
    basename = os.path.basename(input_fa)
    if "_dedup_sorted" in basename:
        id_prefix = basename.split("_dedup_sorted")[0]
    else:
        id_prefix = basename.rsplit(".", 1)[0]
    seq_dict = read_fasta_native(input_fa)
    keep_records = dict()
    for seq_id, seq in seq_dict.items():
        # 给ID加上来源前缀
        prefixed_id = f"{id_prefix}_{seq_id}"
        # 从ID提取count数字，格式 >library_1_Seq_1_reads_10
        count_val = extract_count_from_id(prefixed_id)
        if count_val is None:
            continue
        if count_val >= min_count:
            keep_records[prefixed_id] = seq
    #写入过滤后的fasta
    with open(output_fa, "w", encoding="utf-8") as fw:
        for sid, s in keep_records.items():
            fw.write(f">{sid}\n{s}\n")
    print(f"[预处理过滤] {os.path.basename(input_fa)} 原始:{len(seq_dict)}条,过滤后保留:{len(keep_records)}条(min_count={min_count})")
    return output_fa


def run_cdhit_est(input_fasta, out_prefix, sim, threads=8, n_val=7):
    cmd = [
        CDHIT_PATH,
        "-i", input_fasta,
        "-o", out_prefix,
        "-c", str(sim),
        "-n", str(n_val),
        "-T", str(threads),
        "-M", "0",
        "-d", "0"
    ]
    print(f"执行cd-hit: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    print("====cd-hit stdout====")
    print(result.stdout)
    print("====cd-hit stderr====")
    print(result.stderr)
    if result.returncode != 0:
        print("cd-hit执行返回非0，运行失败！")
        return False
    return True


def parse_cdhit_cluster(clstr_file):
    """解析cd-hit .clst文件，返回家族列表，每个家族包含代表序列，全部成员、总count"""
    families = []
    this_family_members = []
    rep_seqid = None
    with open(clstr_file, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                if len(this_family_members) > 0:
                    families.append({"rep": rep_seqid, "members": this_family_members.copy()})
                this_family_members.clear()
                rep_seqid = None
            elif line != "":
                #解析cd-hit行  0	45nt, >seq1_count100... *
                parts = line.split(", >")
                if len(parts) >= 2:
                    sid_part = parts[-1]
                    seq_id = sid_part.split("...")[0]
                    is_rep = "*" in line
                    this_family_members.append(seq_id)
                    if is_rep:
                        rep_seqid = seq_id
        if len(this_family_members) > 0:
            families.append({"rep": rep_seqid, "members": this_family_members.copy()})
    # 每个家族统计总count
    out_family_list = []
    for fa in families:
        total_cnt = 0
        for mid in fa["members"]:
            c = extract_count_from_id(mid)
            if c is not None:
                total_cnt += c
        out_family_list.append({
            "rep_id": fa["rep"],
            "member_ids": fa["members"],
            "total_count": total_cnt
        })
    #按家族总count降序排序（丰度从高到低）
    out_family_list.sort(key=lambda x: x["total_count"], reverse=True)
    return out_family_list


# ----------------- 步骤2b：家族聚类主流程 -----------------
def run_family_clustering():
    RAW_FASTA_DIR = final_output_dir   # 步骤1的输出目录
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    temp_filter_dir = os.path.join(OUTPUT_DIR, "temp_filter_fa")
    os.makedirs(temp_filter_dir, exist_ok=True)

    # 自动对应 library：
    # 有 Tag -> library_1 ~ library_N
    # 无 Tag -> library_1
    global select_file_names
    if select_file_names is None:
        n_libs = len(second_half_fields) if HAS_TAG else 1
        select_file_names = [
            f"library_{i}_dedup_sorted.fasta"
            for i in range(1, n_libs + 1)
        ]
    # 检查文件是否存在，跳过缺失的文库并给出警告
    existing_file_names = []
    for fname in select_file_names:
        if os.path.exists(os.path.join(RAW_FASTA_DIR, fname)):
            existing_file_names.append(fname)
        else:
            print(f"⚠️ 警告: 未找到 {fname}，该文库将不参与聚类（请检查 tag 数量或步骤1是否成功）")
    if not existing_file_names:
        raise RuntimeError(f"{RAW_FASTA_DIR} 下没有任何 library_xx_dedup_sorted.fasta，无法聚类")
    select_file_names = existing_file_names

    # 预先生成来源前缀列表，供CSV中使用
    id_prefixes = []
    for fname in select_file_names:
        if "_dedup_sorted" in fname:
            id_prefixes.append(fname.split("_dedup_sorted")[0])
        else:
            id_prefixes.append(fname.rsplit(".", 1)[0])

    #1. 循环：对每一个原始fasta做count过滤，生成过滤后的临时fasta
    filtered_fa_list = []
    for fname in select_file_names:
        raw_path = os.path.join(RAW_FASTA_DIR, fname)
        temp_out_fa = os.path.join(temp_filter_dir, f"filtered_{fname}")
        filtered_path = filter_fasta_by_count(raw_path, temp_out_fa, min_count=MIN_RETENTION_COUNT)
        filtered_fa_list.append(filtered_path)

    #2. 将各文库过滤后的fasta合并为1个总fasta
    merged_filtered_fa = os.path.join(OUTPUT_DIR, f"merged_{len(filtered_fa_list)}lib_filtered.fasta")
    with open(merged_filtered_fa, "w", encoding="utf-8") as fw_out:
        for fa in filtered_fa_list:
            with open(fa, "r", encoding="utf-8") as fr:
                fw_out.write(fr.read())
    print(f"\n✅ 各文库合并完成，过滤后总文件：{merged_filtered_fa}")

    #3. cd-hit家族筛选
    cdhit_prefix = os.path.join(OUTPUT_DIR, "cdhit_merged_result")
    ok = run_cdhit_est(merged_filtered_fa, cdhit_prefix, SIM_THRESHOLD, threads=THREADS, n_val=7)
    if not ok:
        raise RuntimeError("cd-hit运行失败")
    clstr_file = cdhit_prefix + ".clstr"
    family_result = parse_cdhit_cluster(clstr_file)

    #====输出文件1：去重全部序列，按家族丰度排序输出fasta====
    all_seq_out = os.path.join(OUTPUT_DIR, "all_seq_sorted_by_family_abund.fasta")
    seq_all_dict = read_fasta_native(merged_filtered_fa)
    with open(all_seq_out, "w", encoding="utf-8") as fout1:
        for fam in family_result:
            for mem_id in fam["member_ids"]:
                seq = seq_all_dict[mem_id]
                fout1.write(f">{mem_id}\n{seq}\n")

    #====输出文件2：丰度最高的各个家族代表序列====
    rep_seq_out = os.path.join(OUTPUT_DIR, "family_representative_top_abund.fasta")
    with open(rep_seq_out, "w", encoding="utf-8") as fout2:
        for idx, fam in enumerate(family_result, start=1):
            rid = fam["rep_id"]
            seq = seq_all_dict[rid]
            fout2.write(f">family_{idx}|{rid}|family_total_count={fam['total_count']}\n{seq}\n")

    #====输出文件3：家族统计表（每成员一行：家族、成员数、来源、序列、counts）====
    stats_out = os.path.join(OUTPUT_DIR, "family_statistics.csv")
    with open(stats_out, "w", encoding="utf-8", newline="") as fcsv:
        writer = csv.writer(fcsv)
        writer.writerow(["Family_Name", "Member_Count", "Source", "Sequence", "Count"])
        for idx, fam in enumerate(family_result, start=1):
            family_name = f"Family_{idx}"
            member_count = len(fam["member_ids"])
            for mid in fam["member_ids"]:
                # 提取来源前缀
                src_prefix = None
                for p in id_prefixes:
                    if mid.startswith(p + "_"):
                        src_prefix = p
                        break
                if src_prefix is None:
                    src_prefix = "unknown"
                # 提取 counts
                c = extract_count_from_id(mid)
                if c is None:
                    c = 0
                seq = seq_all_dict[mid]
                writer.writerow([family_name, member_count, src_prefix, seq, c])

    print("\n====家族聚类任务完成====")
    print(f"输出1【全部去重序列，家族丰度排序】：{all_seq_out}")
    print(f"输出2【家族代表序列，按丰度排序】：{rep_seq_out}")
    print(f"输出3【家族统计表（每成员一行，含成员数）】：{stats_out}")
    return stats_out


# ----------------- 步骤3a：统计表 CSV 转换为家族 txt -----------------
def convert_csv_to_txt(csv_file, output_txt):
    """读取新版 CSV（每成员一行：Family_Name, Member_Count, Source, Sequence, Count），
    按 Family_Name 分组并保持顺序，转换为 motif 代码所需的 txt 格式"""
    os.makedirs(os.path.dirname(output_txt), exist_ok=True)

    with open(csv_file, 'r', encoding='utf-8') as fcsv:
        rows = list(csv.DictReader(fcsv))

    families = []          # [(family_name, [seq1, seq2, ...]), ...]
    cur_name = None
    cur_seqs = []
    for row in rows:
        name = row['Family_Name']
        if cur_name is None:
            cur_name = name
        if name != cur_name:
            families.append((cur_name, cur_seqs))
            cur_name = name
            cur_seqs = []
        seq = row['Sequence'].strip().upper()
        if seq:            # 保留含 N 等 IUPAC 碱基的序列，保证行数与成员数一致
            cur_seqs.append(seq)
    families.append((cur_name, cur_seqs))

    with open(output_txt, 'w', encoding='utf-8') as fout:
        for name, seqs in families:
            fout.write(f"{name}\t{len(seqs)}\n")
            for s in seqs:
                fout.write(f"{s}\n")
            fout.write("\n")  # 家族之间空一行
    print(f"转换完成: {output_txt}，共 {len(families)} 个家族")


# ----------------- 步骤3b：多序列比对（解决 5' 端额外碱基导致的左对齐错位）----------------
def pairwise_align(s1, s2):
    """近似全局比对两条序列。
    返回 (a1, a2, m1, m2)：a1/a2 为等长含 '-' gap 的对齐结果；
    m1[i]/m2[i] 表示 a1[i]/a2[i] 在原序列中的位置，gap 处为 None。"""
    sm = difflib.SequenceMatcher(None, s1, s2, autojunk=False)
    a1, a2, m1, m2 = [], [], [], []
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == 'equal':
            for k in range(i2 - i1):
                a1.append(s1[i1 + k]); a2.append(s2[j1 + k])
                m1.append(i1 + k); m2.append(j1 + k)
        elif tag == 'delete':
            for k in range(i2 - i1):
                a1.append(s1[i1 + k]); a2.append('-')
                m1.append(i1 + k); m2.append(None)
        elif tag == 'insert':
            for k in range(j2 - j1):
                a1.append('-'); a2.append(s2[j1 + k])
                m1.append(None); m2.append(j1 + k)
        elif tag == 'replace':
            l1, l2 = i2 - i1, j2 - j1
            m = max(l1, l2)
            for k in range(m):
                a1.append(s1[i1 + k] if k < l1 else '-')
                a2.append(s2[j1 + k] if k < l2 else '-')
                m1.append(i1 + k if k < l1 else None)
                m2.append(j1 + k if k < l2 else None)
    return ''.join(a1), ''.join(a2), m1, m2


def align_seqs_by_msa(seqs):
    """
    渐进式多序列比对（progressive MSA，列结构实现，正确处理 5'/3' 端任意 indel）。
    原问题：家族序列 5' 端有 0~2 个额外碱基（如 CCACAC... / TCCACAC... / GTCCACAC...），
    直接从头左对齐会使保守核心错位、字母图前端混乱且保守度被拉低。
    改进：以出现次数最多的序列为起点，逐条与共识比对并同步 gap，
    使保守核心真正对齐；5' 端多余碱基变成 gap，不参与核心位点统计。
    仅用 Python 标准库（difflib / collections），无需额外安装。
    返回：等长（含 '-' gap）的对齐序列列表。
    """
    if len(seqs) == 1:
        return list(seqs)
    # 以众数序列为起点，保证起点最"代表"
    cnt = Counter(seqs)
    start = cnt.most_common(1)[0][0]
    n_start = cnt.most_common(1)[0][1]   # 众数序列出现次数（跨文库重复时 >1）
    rest = [s for s in seqs if s != start]
    if not rest:          # 全部是同一序列
        return list(seqs)

    # 列结构：cols[i] = dict {序列索引: 该列字符}；序列无字符则视为 gap
    cols = [{0: ch} for ch in start]
    cons = start
    seq_idx = 1
    for s in rest:
        a_cons, a_s, m_cons, m_s = pairwise_align(cons, s)
        # 以 a_cons 为新的列骨架：m_cons[pos] 给出该列对应 cons 的位置（即 cols 的列号）
        new_cols = []
        for pos in range(len(a_cons)):
            if m_cons[pos] is not None:
                new_cols.append(cols[m_cons[pos]])   # 已有列
            else:
                new_cols.append({})                  # 新列（仅新序列可能有字符）
        # 新序列的字符放入对应列
        for pos in range(len(a_s)):
            if m_s[pos] is not None:
                new_cols[pos][seq_idx] = a_s[pos]
        cols = new_cols
        # 更新共识 = 每列多数投票（含 gap），长度恒等于 len(cols)
        cons_chars = []
        for col in cols:
            if not col:
                cons_chars.append('-')
            else:
                votes = {}
                for ch in col.values():
                    votes[ch] = votes.get(ch, 0) + 1
                cons_chars.append(max(votes, key=votes.get))
        cons = ''.join(cons_chars)
        seq_idx += 1

    # 输出等长（含 '-' gap）的对齐序列
    aligned = []
    for si in range(seq_idx):
        aligned.append(''.join(col.get(si, '-') for col in cols))

    # 【修复】众数序列出现多次时补回 n_start 份：
    # 原实现只保留 1 条众数作对齐起点，其余重复的众数序列被丢弃，
    # 导致 len(aligned) = 1 + len(rest) < len(seqs)，标题 n 比 count 少 n_start-1。
    # 例如同一核心序列在 4 个文库各出现 1 次 → 丢 3 条 → n 比 count 少 3。
    if n_start > 1:
        aligned.extend([aligned[0]] * (n_start - 1))
    return aligned


# ----------------- 步骤3c：绘制 motif 字母图 -----------------
def plot_motif_logos(input_file, csv_file):
    output_dir = os.path.join(OUTPUT_DIR, "motif_logos_info_new")
    png_dir = os.path.join(output_dir, "png")
    pdf_dir = os.path.join(output_dir, "pdf")
    os.makedirs(png_dir, exist_ok=True)
    os.makedirs(pdf_dir, exist_ok=True)

    # 1. 解析家族文件
    families = []
    current = None

    with open(input_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split('\t')
            if parts[0].startswith('Family_'):
                if current is not None:
                    families.append(current)
                current = {
                    'id': parts[0],
                    'count': int(parts[1]) if len(parts) > 1 and parts[1].isdigit() else 0,
                    'seqs': []
                }
            else:
                # 修复：序列行可能没有tab，直接取整行作为序列
                # 接受 IUPAC 碱基（N 等），避免含 N 的序列行被丢弃导致标题 n 与文件 count 不一致
                if current is not None:
                    seq = line.strip().upper()
                    if seq and all(c in 'ATCGNRYKMSWBDHV' for c in seq):
                        current['seqs'].append(seq)
    if current is not None:
        families.append(current)

    print(f"共读取 {len(families)} 个家族")

    # ===== 读取 CSV 中每个家族各成员的 counts（用于筛选"成员少但单序列 counts 高"的家族）=====
    # 适配新版 family_statistics.csv（每成员一行：Family_Name, Source, Sequence, Count）
    fam_member_counts = {}   # Family_X -> [count1, count2, ...]
    if os.path.exists(csv_file):
        with open(csv_file, 'r', encoding='utf-8') as fcsv:
            for row in csv.DictReader(fcsv):
                c = row['Count'].strip()
                if c.isdigit():
                    fam_member_counts.setdefault(row['Family_Name'], []).append(int(c))
    else:
        print(f"[警告] 未找到 CSV: {csv_file}，将只按成员数量>{MEMBER_COUNT_THRESH} 筛选")

    # 按成员数量降序排序，绘制两类家族：
    #   1) 成员数量 > MEMBER_COUNT_THRESH 的家族；
    #   2) 成员数量 < MEMBER_COUNT_THRESH，
    #      但存在某一成员序列 counts > SINGLE_SEQ_COUNT_THRESH 的家族（成员少但单序列丰度高）
    families_sorted = sorted(families, key=lambda x: x['count'], reverse=True)
    target_families = [
        fam for fam in families_sorted
        if fam['count'] > MEMBER_COUNT_THRESH or (
            fam['count'] < MEMBER_COUNT_THRESH and any(c > SINGLE_SEQ_COUNT_THRESH for c in fam_member_counts.get(fam['id'], []))
        )
    ]

    # 分类打印：标注每个家族命中哪种条件
    n_gt50 = sum(1 for fam in target_families if fam['count'] > MEMBER_COUNT_THRESH)
    n_lt50_high = len(target_families) - n_gt50
    print(f"将绘制 {len(target_families)} 个家族（成员数>{MEMBER_COUNT_THRESH}: {n_gt50} 个；"
          f"成员数<{MEMBER_COUNT_THRESH} 且单序列counts>{SINGLE_SEQ_COUNT_THRESH}: {n_lt50_high} 个）:")
    for fam in target_families:
        if fam['count'] > MEMBER_COUNT_THRESH:
            tag = f"成员数>{MEMBER_COUNT_THRESH}"
        else:
            maxc = max(fam_member_counts.get(fam['id'], [0]) or [0])
            tag = f"成员数<{MEMBER_COUNT_THRESH}且单序列counts>{SINGLE_SEQ_COUNT_THRESH} (max={maxc})"
        print(f"  {fam['id']}: {fam['count']} 条序列  [{tag}]")

    for fam in target_families:
        fam_id = fam['id']
        seqs = fam['seqs']

        if len(seqs) < 2:
            print(f"跳过 {fam_id}：只有 {len(seqs)} 条序列")
            continue

        # ===== 先做多序列比对再统计（解决 5' 端额外碱基导致的左对齐错位）=====
        # 直接左对齐会使保守核心错位、字母图前端混乱、整体保守度被拉低；截断还会丢掉尾部保守位点。
        seqs = align_seqs_by_msa(seqs)
        seq_len = len(seqs[0])

        if seq_len == 0:
            print(f"跳过 {fam_id}：序列长度为 0")
            continue

        # ========== 统计概率矩阵（跳过 gap，按该列 ATCG 总数归一）==========
        counts = {base: [0] * seq_len for base in 'ATCG'}
        col_cov = [0] * seq_len        # 每列非 gap 碱基数（覆盖度）
        for seq in seqs:
            for i, base in enumerate(seq):
                if base in counts:
                    counts[base][i] += 1
                    col_cov[i] += 1

        n_seqs = len(seqs)
        # 丢弃覆盖度过低的列（5' 端只有少数序列才有的前缀位），避免小样本假信号
        keep_cols = [i for i in range(seq_len) if col_cov[i] >= max(int(MIN_COL_COVERAGE * n_seqs), 3)]

        if not keep_cols:
            print(f"跳过 {fam_id}：所有列覆盖度过低")
            continue

        df_prob = pd.DataFrame({
            base: [counts[base][i] / max(col_cov[i], 1) for i in keep_cols]
            for base in 'ATCG'
        })
        seq_len = len(keep_cols)

        # ========== 概率 → 信息含量 (bits) ==========
        entropy = np.zeros(seq_len)
        for i in range(seq_len):
            for base in 'ATCG':
                p = df_prob.iloc[i][base]
                if p > 0:
                    entropy[i] -= p * np.log2(p)

        info_content = 2.0 - entropy  # 最大信息含量 = log2(4) = 2 bits

        # 每个碱基的标准 sequence logo 高度 = p_i * IC_i
        df_logo = pd.DataFrame({
            base: [df_prob.iloc[i][base] * info_content[i] for i in range(seq_len)]
            for base in 'ATCG'
        })

        # ========== 绘图 ==========
        fig_width = max(seq_len * 0.4, 5)
        fig, ax = plt.subplots(figsize=(fig_width, 3.8))

        logo = logomaker.Logo(
            df_logo,
            ax=ax,
            color_scheme=DNA_COLORS,
            font_name='DejaVu Sans',
            vpad=0.05,
            width=0.9
        )

        # 字母加粗 + 细黑边
        logo.style_glyphs(fontweight='bold', stroke_width=0.5, stroke_color='black')

        # 美化边框
        logo.style_spines(visible=False)
        logo.style_spines(spines=['left', 'bottom'], visible=True, linewidth=1.2)
        ax.spines['left'].set_bounds(0, 2.0)
        ax.spines['bottom'].set_bounds(0, seq_len - 1)

        # 标题与轴标签
        ax.set_title(f"{fam_id}  (n={n_seqs} sequences)", fontsize=14, fontweight='bold', pad=12)
        ax.set_xlabel("Position (bp)", fontsize=12)
        ax.set_ylabel("Information content (bits)", fontsize=12)

        # 纵轴固定 0~2
        ax.set_ylim(0, 2.0)
        ax.set_yticks([0, 0.5, 1.0, 1.5, 2.0])
        ax.set_yticklabels(['0', '0.5', '1.0', '1.5', '2.0'], fontsize=10)
        ax.axhline(y=2.0, color='gray', linestyle='--', linewidth=0.8, alpha=0.5)

        # X轴自适应
        if seq_len <= 30:
            step = 1
        elif seq_len <= 60:
            step = 2
        elif seq_len <= 100:
            step = 5
        else:
            step = 10

        ax.set_xticks(range(0, seq_len, step))
        ax.set_xticklabels(range(1, seq_len + 1, step), rotation=90 if seq_len > 40 else 0, fontsize=9)
        ax.tick_params(axis='both', which='major', labelsize=10)

        plt.tight_layout()

        # 保存 PNG（png/ 文件夹）和 PDF（pdf/ 文件夹）
        out_png = os.path.join(png_dir, f"{fam_id}_motif.png")
        out_pdf = os.path.join(pdf_dir, f"{fam_id}_motif.pdf")
        plt.savefig(out_png, dpi=300, bbox_inches='tight', facecolor='white')
        plt.savefig(out_pdf, bbox_inches='tight', facecolor='white')
        plt.close()

        print(f"✓ 已保存: {out_png}  +  {out_pdf}")

    print(f"\n全部完成！PNG 保存在: {png_dir}\nPDF 保存在: {pdf_dir}")


# ================================ 主流程 ================================
if __name__ == "__main__":
    print("🚀 开始 SELEX 文库序列完整分析流程...")
    print("=" * 72)
    print(f"模式：{'有 Tag' if HAS_TAG else '无 Tag'}")
    print(f"random_length：{random_length}（完全由用户配置，与 Tag 模式无关）")
    if HAS_TAG:
        print(f"Tag 数量：{len(second_half_fields)}")
        print(f"Tag：{second_half_fields}")
    else:
        print("Tag：未填写，所有序列归入 library_1")
    print("=" * 72)

    # 步骤1：去接头 → FLASH拼接 → （可选）去PCR引物 → 提取随机区 → 文库归类
    process_pipeline()

    # 步骤2：count过滤 → 合并 → cd-hit家族聚类 → 输出 fasta + CSV
    stats_csv = run_family_clustering()

    # 步骤3a：CSV 转 txt
    fam_txt = os.path.join(OUTPUT_DIR, "dna_families_with_seq.txt")
    convert_csv_to_txt(stats_csv, fam_txt)

    # 步骤3b：多序列比对 + 绘制 motif 字母图
plot_motif_logos(fam_txt, stats_csv)