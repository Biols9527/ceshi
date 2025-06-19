#!/usr/bin/env python3
import sys
import os
import getopt
sys.path.append('/dellfsqd2/ST_OCEAN/USER/lishuo1/01_software/miniconda3/lib/python3.7/site-packages/')
import subprocess
import argparse
import logging
import traceback
from concurrent.futures import ThreadPoolExecutor
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from functools import lru_cache
import re
import time
import shutil
import uuid
import json
import itertools
import hashlib
import multiprocessing
# 提前导入可能需要的所有模块
from Bio.Blast import NCBIXML
import numpy as np
from functools import lru_cache

# 设置日志记录，带有时间戳
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.FileHandler("msa_pipeline.log"),
        logging.StreamHandler()
    ]
)
logger = logging.getLogger(__name__)

# 定义4倍简并密码子（第三位可以是任何核苷酸）
FOURFOLD_CODONS = {
    'GC': True,  # 丙氨酸 (GCN)
    'CG': True,  # 精氨酸 (CGN)
    'GG': True,  # 甘氨酸 (GGN)
    'CT': True,  # 亮氨酸 (CTN)
    'CC': True,  # 脯氨酸 (CCN)
    'TC': True,  # 丝氨酸 (TCN)
    'AC': True,  # 苏氨酸 (ACN)
    'GT': True   # 缬氨酸 (GTN)
}

# 定义遗传密码表，将密码子映射到氨基酸（包括模糊密码子）
GENETIC_CODE = {
    'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
    'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
    'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
    'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
    'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
    'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
    'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
    'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
    'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
    'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
    'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
    'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
    'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
}

# 最小报告重复次数 - 提高阈值减少误报
MIN_REPEAT_COUNT_TO_REPORT = 70  # 提高阈值，减少CDS中正常重复的警告

# 性能优化：BLAST数据库缓存
BLAST_DB_CACHE = {}
BLAST_DB_CACHE_SIZE = 5  # 最大缓存数据库数量

# 性能优化：缓存函数计算结果
@lru_cache(maxsize=1024)
def safe_translate_codon_cached(codon):
    """
    安全地将密码子翻译为氨基酸的缓存版本
    
    Args:
        codon: 三字母核苷酸密码子
        
    Returns:
        单字母氨基酸或'X'表示未知
    """
    return safe_translate_codon(codon)

def find_cds_files(directory='.'):
    """查找给定目录中的所有.cds文件"""
    try:
        cds_files = [f for f in os.listdir(directory) if f.endswith('.cds')]
        logger.info(f"找到 {len(cds_files)} 个 .cds 文件在 {directory} 中")
        return sorted(cds_files)  # 排序，保证处理顺序一致
    except Exception as e:
        logger.error(f"在 {directory} 中查找 CDS 文件时出错: {str(e)}")
        return []

def extract_species_name(header):
    """
    根据FASTA头部提取物种名称，支持">ID SPECIES_NAME"格式
    
    Args:
        header: FASTA头部字符串
    
    Returns:
        提取的物种名称
    """
    # 移除'>'前缀如果存在
    if header.startswith('>'):
        header = header[1:]
    
    # 按空格分割，第一部分是ID，后面的是物种名
    parts = header.split(' ', 1)
    
    if len(parts) > 1:
        # 返回空格后的部分作为物种名
        return parts[1].strip()
    else:
        # 如果没有空格，则使用整个标识符作为物种名
        return header.strip()

def check_coding_sequence(sequence):
    """
    检查序列是否看起来像有效的编码序列
    
    Returns:
        tuple: (is_valid, message)
    """
    if not sequence:
        return False, "空序列"
        
    if len(sequence) % 3 != 0:
        return False, f"序列长度 {len(sequence)} 不是 3 的倍数"
    
    # 检查中间的终止密码子
    for i in range(0, len(sequence) - 3, 3):  # 排除最后一个密码子
        codon = sequence[i:i+3].upper()
        if codon in ['TAA', 'TAG', 'TGA'] and i < len(sequence)-5:
            return False, f"序列中间发现终止密码子 {codon} 在位置 {i}"
    
    # 检查无效字符
    valid_chars = set('ACGTRYKMSWBDHVN-')  # DNA + 模糊代码 + 缺口
    invalid_chars = set(sequence.upper()) - valid_chars
    if invalid_chars:
        return False, f"序列包含无效字符: {', '.join(invalid_chars)}"
    
    return True, "有效的编码序列"

def analyze_codon_structure(sequence):
    """
    分析序列的密码子结构，识别可能的起始/终止位点和开放阅读框
    
    Args:
        sequence: 核苷酸序列
    
    Returns:
        包含分析结果的字典
    """
    # 性能优化：使用numpy数组加速处理
    # 去除缺口以进行更准确的ORF分析
    clean_seq = sequence.upper().replace('-', '').replace('N', '')
    
    # 分析3个可能的阅读框
    frames = []
    frame_orfs = []
    
    for frame in range(3):
        codons = []
        for i in range(frame, len(clean_seq)-2, 3):
            codons.append(clean_seq[i:i+3])
        frames.append(codons)
        
        # 查找这个阅读框中的所有ORFs
        orfs = []
        current_orf = []
        
        for i, codon in enumerate(codons):
            if codon in ['TAA', 'TAG', 'TGA']:
                if current_orf:
                    orfs.append({
                        'start': frame + 3*frames[frame].index(current_orf[0]),
                        'end': frame + 3*(frames[frame].index(current_orf[0])+len(current_orf)),
                        'codons': current_orf
                    })
                    current_orf = []
            else:
                current_orf.append(codon)
        
        # 添加最后一个ORF(如果没有终止密码子)
        if current_orf:
            orfs.append({
                'start': frame + 3*frames[frame].index(current_orf[0]),
                'end': frame + 3*(frames[frame].index(current_orf[0])+len(current_orf)),
                'codons': current_orf
            })
        
        frame_orfs.append(orfs)
    
    # 查找起始密码子
    start_positions = []
    for frame_idx, frame_codons in enumerate(frames):
        for i, codon in enumerate(frame_codons):
            if codon == 'ATG':
                # 记录可能的起始位置
                start_positions.append({
                    'frame': frame_idx,
                    'position': frame_idx + i*3,
                    'codon_idx': i
                })
    
    # 找最长的ORF
    longest_orf = None
    best_frame = 0
    
    for frame_idx, orfs in enumerate(frame_orfs):
        for orf in orfs:
            if not longest_orf or len(orf['codons']) > len(longest_orf['codons']):
                longest_orf = orf
                best_frame = frame_idx
    
    # 检查最长ORF是否以起始密码子开始
    has_start = False
    if longest_orf and longest_orf['codons']:
        has_start = (longest_orf['codons'][0] == 'ATG')
    
    result = {
        'best_frame': best_frame,
        'longest_orf': longest_orf,
        'has_start': has_start,
        'has_stop': False,  # 最长ORF可能没有终止子，会在后续逻辑中处理
        'start_positions': start_positions,
        'orf_length': len(longest_orf['codons'])*3 if longest_orf else 0
    }
    
    return result

def assess_sequence_quality(sequence, seq_id="Unknown"):
    """
    评估序列质量，检测潜在问题
    仅报告非常极端的重复情况，忽略CDS中正常的重复
    
    Args:
        sequence: 核苷酸序列
        seq_id: 序列ID (用于日志)
    
    Returns:
        问题列表
    """
    issues = []
    
    # 清理序列以便更准确的分析
    clean_seq = sequence.upper().replace('-', '')
    
    # 检查序列长度
    if len(clean_seq) < 30:  # 小于10个密码子
        issues.append("序列过短")
    
    # 检查GC含量
    gc_count = clean_seq.count('G') + clean_seq.count('C')
    gc_content = gc_count / len(clean_seq) if clean_seq else 0
    if gc_content < 0.25 or gc_content > 0.75:
        issues.append(f"异常GC含量: {gc_content:.2f}")
    
    # 检查未知碱基比例
    n_count = clean_seq.count('N')
    n_ratio = n_count / len(clean_seq) if clean_seq else 0
    if n_ratio > 0.1:
        issues.append(f"高比例未知碱基: {n_ratio:.2f}")
    
    # 检查终止密码子
    if len(clean_seq) % 3 == 0:  # 只有当序列是3的倍数时才检查密码子
        for i in range(0, len(clean_seq)-5, 3):  # 排除最后一个密码子
            codon = clean_seq[i:i+3]
            if codon in ['TAA', 'TAG', 'TGA'] and i < len(clean_seq)-5:
                issues.append(f"序列中部有终止密码子，位置 {i}")
    
    # 性能优化：只检查7和10碱基长度的重复，这不是3的倍数，可能暗示frame shift问题
    logger.debug(f"分析序列 {seq_id} 的重复模式 (长度: {len(clean_seq)})")

    repeat_checks = [7, 10]  # 非3倍数长度，可能表示框架问题
    significant_repeats = []
    
    # 性能优化：只在序列长度足够时才检查重复
    if len(clean_seq) > 100:
        for repeat_len in repeat_checks:
            # 优化：使用字典直接计数而不是构建中间列表
            repeat_counts = defaultdict(int)
            
            for i in range(len(clean_seq) - repeat_len + 1):
                fragment = clean_seq[i:i+repeat_len]
                if 'N' not in fragment:  # 忽略含N的片段
                    repeat_counts[fragment] += 1
                    
            # 计算总重复数，只统计出现超过1次的片段
            total_repeats = sum(1 for count in repeat_counts.values() if count > 1)
            
            # 使用更高的阈值，减少误报
            threshold = MIN_REPEAT_COUNT_TO_REPORT
            if repeat_len == 10:
                threshold = MIN_REPEAT_COUNT_TO_REPORT // 2  # 较长重复的阈值降低
                
            # 只有当重复数大幅超过阈值时才报告，这可能表示真正的问题
            if total_repeats >= threshold:
                significant_repeats.append(f"发现 {total_repeats} 个长度为 {repeat_len} 的重复 (非3倍数，可能影响密码子框架)")
                logger.debug(f"序列 {seq_id} 中{repeat_len}碱基重复统计: 总计 {total_repeats} 个")
    
    # 只有存在显著重复时才添加到问题列表
    if significant_repeats:
        issues.extend(significant_repeats)
    
    return issues

def handle_non_triplet_sequence(sequence, seq_id, position='end', context=None):
    """
    更智能地处理非3倍数的序列
    
    Args:
        sequence: 原始序列
        seq_id: 序列ID用于日志
        position: 填充位置 'start', 'end', 'smart'
        context: 提供额外上下文信息的词典
    
    Returns:
        处理后的序列
    """
    remainder = len(sequence) % 3
    if remainder == 0:
        return sequence
    
    padding_needed = 3 - remainder
    
    # 根据序列特征决定填充策略
    if position == 'smart' and context is not None:
        # 检查序列是否有起始/终止密码子
        has_start = context.get('has_start', False) or sequence[:3].upper() == 'ATG'
        has_stop = context.get('has_stop', False) or sequence[-3:].upper() in ['TAA', 'TAG', 'TGA']
        
        # 基于密码子情况决定填充位置
        if has_start and not has_stop:
            # 可能是3'端不完整
            padding = 'N' * padding_needed
            padded_seq = sequence + padding
            logger.debug(f"序列 {seq_id} 有起始密码子但无终止密码子，在3'端添加 {padding_needed} 个N")
        elif not has_start and has_stop:
            # 可能是5'端不完整
            padding = 'N' * padding_needed
            padded_seq = padding + sequence
            logger.debug(f"序列 {seq_id} 无起始密码子但有终止密码子，在5'端添加 {padding_needed} 个N")
        else:
            # 默认在末尾添加，但记录更详细的信息
            padding = 'N' * padding_needed
            padded_seq = sequence + padding
            logger.debug(f"序列 {seq_id} 不含完整密码子信息，在3'端添加 {padding_needed} 个N")
    elif position == 'start':
        padding = 'N' * padding_needed
        padded_seq = padding + sequence
        logger.debug(f"序列 {seq_id} 在5'端添加 {padding_needed} 个N")
    else:  # 'end' 或默认
        padding = 'N' * padding_needed
        padded_seq = sequence + padding
        logger.debug(f"序列 {seq_id} 在3'端添加 {padding_needed} 个N")
    
    return padded_seq

def cleanup_temp_files(temp_fasta, temp_output, db_prefix):
    """清理临时BLAST文件，但保留缓存的数据库"""
    try:
        # 删除临时FASTA和输出文件
        for file in [temp_fasta, temp_output]:
            if file and os.path.exists(file):
                os.remove(file)
                logger.debug(f"已删除临时文件: {file}")
        
        # 不删除可能被缓存的数据库文件
        if db_prefix and db_prefix not in BLAST_DB_CACHE:
            for ext in [".nhr", ".nin", ".nsq", ".fasta"]:
                db_file = db_prefix + ext
                if os.path.exists(db_file):
                    os.remove(db_file)
                    logger.debug(f"已删除BLAST数据库文件: {db_file}")
    except Exception as e:
        logger.debug(f"清理临时文件时出错: {str(e)}")

def hash_file_content(file_path):
    """
    计算文件内容的哈希值，用于数据库缓存
    """
    try:
        hasher = hashlib.md5()
        with open(file_path, 'rb') as f:
            buf = f.read(65536)  # 64kb chunks
            while buf:
                hasher.update(buf)
                buf = f.read(65536)
        return hasher.hexdigest()
    except Exception as e:
        logger.debug(f"计算文件哈希时出错: {str(e)}")
        return None

def create_or_get_blast_db(file_path, temp_dir, blast_path):
    """
    创建或从缓存获取BLAST数据库
    
    Args:
        file_path: CDS文件路径
        temp_dir: 临时目录
        blast_path: BLAST可执行文件路径
    
    Returns:
        临时数据库路径或None
    """
    # 计算文件内容哈希作为缓存键
    file_hash = hash_file_content(file_path)
    if not file_hash:
        return create_temp_blast_db(file_path, temp_dir, blast_path)
    
    # 检查缓存中是否已存在数据库
    if file_hash in BLAST_DB_CACHE:
        db_prefix = BLAST_DB_CACHE[file_hash]
        # 验证数据库文件是否仍然存在
        if os.path.exists(db_prefix + ".nsq"):
            logger.info(f"使用缓存的BLAST数据库: {db_prefix}")
            return db_prefix
        else:
            # 缓存项有效但文件已删除，从缓存中移除
            del BLAST_DB_CACHE[file_hash]
    
    # 创建新数据库
    db_prefix = create_temp_blast_db(file_path, temp_dir, blast_path)
    if db_prefix:
        # 管理缓存大小
        if len(BLAST_DB_CACHE) >= BLAST_DB_CACHE_SIZE:
            # 删除最旧的数据库
            oldest_key = next(iter(BLAST_DB_CACHE))
            oldest_db = BLAST_DB_CACHE[oldest_key]
            del BLAST_DB_CACHE[oldest_key]
            # 清理关联文件
            for ext in [".nhr", ".nin", ".nsq", ".fasta"]:
                old_file = oldest_db + ext
                if os.path.exists(old_file):
                    os.remove(old_file)
        
        # 添加新数据库到缓存
        BLAST_DB_CACHE[file_hash] = db_prefix
    
    return db_prefix

def create_temp_blast_db(file_path, temp_dir, blast_path):
    """
    从输入文件自动创建临时BLAST数据库
    
    Args:
        file_path: CDS文件路径
        temp_dir: 临时目录
        blast_path: BLAST可执行文件路径
    
    Returns:
        临时数据库路径或None
    """
    try:
        # 确保临时目录存在
        os.makedirs(temp_dir, exist_ok=True)
        
        logger.info(f"开始为 {os.path.basename(file_path)} 创建临时BLAST数据库")
        
        # 读取所有有效序列（3的倍数）- 性能优化：批量读取
        valid_seqs = []
        total_seqs = 0
        
        # 使用生成器而不是一次读取所有
        for record in SeqIO.parse(file_path, "fasta"):
            total_seqs += 1
            sequence = str(record.seq)
            if len(sequence) % 3 == 0:
                valid_seqs.append(record)
                # 当收集足够的序列后，停止读取
                if len(valid_seqs) >= 100:  # 限制使用的序列数量
                    break
        
        # 如果没有有效序列，返回None
        if not valid_seqs:
            logger.warning(f"文件 {file_path} 中没有长度为3倍数的有效序列, 无法创建BLAST数据库")
            return None
        
        logger.info(f"从 {file_path} 中找到 {len(valid_seqs)}/{total_seqs} 个有效序列用于创建BLAST数据库")
        
        # 创建临时FASTA文件
        db_prefix = os.path.join(temp_dir, f"temp_blastdb_{uuid.uuid4()}")
        temp_fasta = f"{db_prefix}.fasta"
        
        with open(temp_fasta, 'w') as f:
            for record in valid_seqs:
                f.write(f">{record.id}\n{record.seq}\n")
        
        # 构建makeblastdb路径
        blast_dir = os.path.dirname(blast_path)
        makeblastdb_path = os.path.join(blast_dir, "makeblastdb")
        
        # 尝试几种可能的路径
        makeblastdb_candidates = [
            makeblastdb_path,
            blast_path.replace("blastn", "makeblastdb"),
            os.path.join(blast_dir, "makeblastdb.exe"),  # Windows
            blast_path.replace("blastn.exe", "makeblastdb.exe")  # Windows
        ]
        
        makeblastdb_path = None
        for candidate in makeblastdb_candidates:
            if os.path.exists(candidate):
                makeblastdb_path = candidate
                break
                
        if not makeblastdb_path:
            logger.error(f"无法找到makeblastdb可执行文件，尝试了: {', '.join(makeblastdb_candidates)}")
            return None
            
        logger.debug(f"使用makeblastdb路径: {makeblastdb_path}")
        
        # 运行makeblastdb
        cmd = [
            makeblastdb_path,
            "-in", temp_fasta,
            "-dbtype", "nucl",
            "-out", db_prefix
        ]
        
        logger.debug(f"创建BLAST数据库命令: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.error(f"创建BLAST数据库失败: {result.stderr}")
            return None
        
        logger.info(f"成功创建临时BLAST数据库: {db_prefix}")
        return db_prefix
    
    except Exception as e:
        logger.error(f"创建临时BLAST数据库时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def handle_non_triplet_with_blast(sequence, seq_id, file_path, blast_path, temp_dir):
    """
    使用BLAST处理非3倍数序列，修改逻辑以明确仅在必要时调用
    性能优化：共享和缓存BLAST数据库
    
    Args:
        sequence: 原始序列
        seq_id: 序列ID
        file_path: 原始CDS文件路径
        blast_path: BLAST可执行文件路径
        temp_dir: 临时文件目录
    
    Returns:
        处理后的序列
    """
    # 只对非3倍数序列启动BLAST处理
    if len(sequence) % 3 == 0:
        return sequence
    
    logger.info(f"启动BLAST处理非3倍数序列: {seq_id}, 长度 {len(sequence)}")
        
    if not blast_path or not temp_dir or not file_path:
        logger.warning(f"未提供完整BLAST参数，使用默认填充方法处理序列 {seq_id}")
        return handle_non_triplet_sequence(sequence, seq_id, position='smart')
    
    try:
        # 获取或创建BLAST数据库
        db_prefix = create_or_get_blast_db(file_path, temp_dir, blast_path)
        
        if not db_prefix:
            logger.warning(f"无法创建临时BLAST数据库，使用默认填充方法处理序列 {seq_id}")
            return handle_non_triplet_sequence(sequence, seq_id, position='smart')
        
        # 创建临时序列文件
        temp_fasta = os.path.join(temp_dir, f"{seq_id}_{uuid.uuid4()}.fasta")
        
        with open(temp_fasta, 'w') as f:
            f.write(f">{seq_id}\n{sequence}\n")
        
        # 运行BLASTN
        temp_output = os.path.join(temp_dir, f"{seq_id}_{uuid.uuid4()}.xml")
        
        cmd = [
            blast_path,
            "-db", db_prefix,
            "-query", temp_fasta,
            "-outfmt", "5",  # XML格式输出
            "-out", temp_output,
            "-task", "blastn",  # 使用普通blastn而不是megablast以提高敏感度
            "-word_size", "7",  # 使用较小的word_size以便匹配更短的片段
            "-evalue", "0.01"   # 较严格的E值阈值
        ]
        
        logger.debug(f"BLASTN命令: {' '.join(cmd)}")
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            logger.warning(f"BLASTN运行失败: {result.stderr}")
            # 清理临时文件
            cleanup_temp_files(temp_fasta, temp_output, None)  # 不删除数据库
            return handle_non_triplet_sequence(sequence, seq_id, position='smart')
        
        # 解析BLASTN结果
        best_hit = None
        with open(temp_output) as f:
            blast_records = NCBIXML.parse(f)
            for blast_record in blast_records:
                if blast_record.alignments:
                    best_hit = blast_record.alignments[0]
                    break
        
        if not best_hit:
            logger.warning(f"序列 {seq_id} 未找到BLAST匹配结果")
            # 清理临时文件
            cleanup_temp_files(temp_fasta, temp_output, None)
            return handle_non_triplet_sequence(sequence, seq_id, position='smart')
        
        # 分析比对结果
        hsp = best_hit.hsps[0]  # 最佳高分片段对
        
        # 检查框架和填充位置
        query_start = hsp.query_start  # 查询序列起始位置
        query_end = hsp.query_end  # 查询序列结束位置
        sbjct_start = hsp.sbjct_start  # 匹配序列起始位置
        sbjct_end = hsp.sbjct_end  # 匹配序列结束位置
        
        # 根据比对结果确定填充策略
        remainder = len(sequence) % 3
        padding_needed = 3 - remainder
        
        if query_start > 3:  # 序列前端有超过3个碱基的未比对部分
            padding = 'N' * padding_needed
            padded_seq = padding + sequence
            logger.info(f"序列 {seq_id} 根据BLAST结果在5'端添加 {padding_needed} 个N (起始位置: {query_start})")
        elif len(sequence) - query_end > 3:  # 序列后端有超过3个碱基的未比对部分
            padding = 'N' * padding_needed
            padded_seq = sequence + padding
            logger.info(f"序列 {seq_id} 根据BLAST结果在3'端添加 {padding_needed} 个N (结束位置: {query_end}/{len(sequence)})")
        else:
            # 根据主要比对对象的阅读框架推断
            frame_offset = (sbjct_start - 1) % 3
            current_offset = (query_start - 1) % 3
            
            if frame_offset != current_offset:
                # 需要调整框架使其与参考一致
                if current_offset < frame_offset:
                    padding = 'N' * (frame_offset - current_offset)
                    padded_seq = padding + sequence
                    logger.info(f"序列 {seq_id} 根据参考阅读框在5'端添加 {len(padding)} 个N (框架调整)")
                else:
                    padding = 'N' * (3 - (current_offset - frame_offset))
                    padded_seq = padding + sequence
                    logger.info(f"序列 {seq_id} 根据参考阅读框在5'端添加 {len(padding)} 个N (框架调整)")
            else:
                # 框架一致，根据序列内容分析添加填充
                # 检查是否有明确的起始或终止密码子
                has_start = sequence[:3].upper() == 'ATG'
                has_stop = sequence[-3:].upper() in ['TAA', 'TAG', 'TGA']
                
                if has_start and not has_stop:
                    # 有起始密码子但没有终止密码子，可能是3'端缺失
                    padding = 'N' * padding_needed
                    padded_seq = sequence + padding
                    logger.info(f"序列 {seq_id} 有起始密码子但无终止密码子，在3'端添加 {padding_needed} 个N")
                elif not has_start and has_stop:
                    # 没有起始密码子但有终止密码子，可能是5'端缺失
                    padding = 'N' * padding_needed
                    padded_seq = padding + sequence
                    logger.info(f"序列 {seq_id} 无起始密码子但有终止密码子，在5'端添加 {padding_needed} 个N")
                else:
                    # 默认在序列末尾添加填充
                    padding = 'N' * padding_needed
                    padded_seq = sequence + padding
                    logger.info(f"序列 {seq_id} 框架与参考一致，在3'端添加 {padding_needed} 个N (默认策略)")
        
        # 清理临时文件
        cleanup_temp_files(temp_fasta, temp_output, None)  # 不删除数据库
        
        return padded_seq
        
    except Exception as e:
        logger.error(f"使用BLAST处理序列 {seq_id} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return handle_non_triplet_sequence(sequence, seq_id, position='smart')

def preprocess_cds_sequence(sequence, seq_id, file_path=None, blast_params=None):
    """
    序列预处理流水线，包括质量检查、结构分析和必要的修复
    修改逻辑，只对非3倍数序列使用BLAST
    
    Args:
        sequence: 原始序列
        seq_id: 序列ID
        file_path: 原始CDS文件路径（用于创建BLAST数据库）
        blast_params: BLAST参数字典
    
    Returns:
        处理后的序列和结构信息
    """
    # 性能优化：仅在调试模式下执行详细的质量检查
    if logger.getEffectiveLevel() <= logging.DEBUG:
        issues = assess_sequence_quality(sequence, seq_id)
        if issues:
            logger.warning(f"序列 {seq_id} 存在质量问题: {', '.join(issues)}")
    
    # 分析密码子结构
    structure = analyze_codon_structure(sequence)
    
    # 如果需要，调整序列框架
    if structure['best_frame'] != 0:
        logger.info(f"序列 {seq_id} 最佳阅读框为 {structure['best_frame']}，调整中")
        sequence = sequence[structure['best_frame']:]
    
    # 处理非3倍数长度
    if len(sequence) % 3 != 0:
        # 只有非3倍数序列才需要处理
        logger.info(f"发现非3倍数序列: {seq_id}, 长度 {len(sequence)}")
        
        if blast_params and blast_params.get('use_blast', False) and file_path:
            # 使用BLAST处理
            logger.info(f"使用BLAST处理非3倍数序列 {seq_id}")
            sequence = handle_non_triplet_with_blast(
                sequence, 
                seq_id,
                file_path,
                blast_params.get('blast_path'),
                blast_params.get('temp_dir')
            )
        else:
            # 使用基于序列特征的处理
            context = {
                'has_start': structure.get('has_start', False),
                'has_stop': structure.get('has_stop', False),
                'best_frame': structure.get('best_frame', 0)
            }
            sequence = handle_non_triplet_sequence(sequence, seq_id, position='smart', context=context)
    
    return sequence, structure

def parse_fasta(file_path, duplicate_strategy='longest', blast_params=None):
    """
    从FASTA文件解析序列，使用指定的重复处理策略
    修复物种名称提取和BLAST处理逻辑
    性能优化：批量处理序列
    
    Args:
        file_path: FASTA文件路径
        duplicate_strategy: 重复处理策略
        blast_params: BLAST参数
    
    Returns:
        物种到序列的映射字典
    """
    logger.debug(f"解析 FASTA 文件: {file_path}, 策略: {duplicate_strategy}")
    
    # 如果策略是'alignment_quality'，使用专门的函数处理
    if duplicate_strategy == 'alignment_quality':
        return parse_fasta_with_duplicates(file_path, blast_params=blast_params)
    
    species_seqs = {}
    species_counts = defaultdict(int)
    duplicate_found = False
    original_ids = {}  # 存储原始ID作为参考
    
    try:
        # 性能优化：使用生成器而不是一次性加载所有记录
        record_count = 0
        batch_size = 100  # 每批处理的序列数
        batch = []
        
        for record in SeqIO.parse(file_path, "fasta"):
            batch.append(record)
            record_count += 1
            
            # 批量处理
            if len(batch) >= batch_size:
                process_fasta_batch(batch, species_seqs, species_counts, original_ids, duplicate_strategy, file_path, blast_params)
                batch = []
                logger.debug(f"已处理 {record_count} 条序列")
        
        # 处理最后一批
        if batch:
            process_fasta_batch(batch, species_seqs, species_counts, original_ids, duplicate_strategy, file_path, blast_params)
            
        logger.debug(f"在 {file_path} 中找到并处理了 {record_count} 个序列")
        
        # 检查是否有重复
        duplicate_found = any(count > 1 for count in species_counts.values())
        if duplicate_found:
            logger.info(f"在 {file_path} 中发现重复物种: {dict(species_counts)}")
        
        return species_seqs
    except Exception as e:
        logger.error(f"解析 FASTA 文件 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def process_fasta_batch(batch, species_seqs, species_counts, original_ids, duplicate_strategy, file_path, blast_params):
    """
    批量处理FASTA记录
    
    Args:
        batch: 要处理的记录批次
        species_seqs: 物种到序列的映射
        species_counts: 物种计数字典
        original_ids: 原始ID映射
        duplicate_strategy: 重复处理策略
        file_path: 原始文件路径
        blast_params: BLAST参数
    """
    for record in batch:
        # 从头部提取物种名
        species = extract_species_name(record.description)
        sequence = str(record.seq)
        
        # 预处理序列
        sequence, _ = preprocess_cds_sequence(sequence, record.id, file_path, blast_params)
        
        # 检查序列是否是有效的编码序列
        is_valid, message = check_coding_sequence(sequence)
        if not is_valid:
            logger.warning(f"序列 {record.id} ({species}) 可能不是有效的编码序列: {message}")
        
        # 检查这个物种是否已存在
        if species in species_seqs:
            species_counts[species] += 1
            
            if duplicate_strategy == 'longest':
                # 保留最长的序列
                if len(sequence) > len(species_seqs[species]):
                    logger.debug(f"用更长的序列替换 {species} 的较短序列 (长度: {len(species_seqs[species])} → {len(sequence)})")
                    species_seqs[species] = sequence
                    original_ids[species] = record.id
            
            elif duplicate_strategy == 'first':
                # 保留第一个序列，忽略这个
                logger.debug(f"忽略 {species} 的重复序列")
                continue
                
            elif duplicate_strategy == 'rename':
                # 通过添加后缀重命名重复物种
                new_species = f"{species}_{species_counts[species]}"
                logger.debug(f"将重复的物种从 {species} 重命名为 {new_species}")
                species_seqs[new_species] = sequence
                original_ids[new_species] = record.id
        else:
            # 第一次见到这个物种
            species_seqs[species] = sequence
            species_counts[species] = 1
            original_ids[species] = record.id

def parse_fasta_with_duplicates(file_path, blast_params=None):
    """
    解析FASTA文件，保留所有重复的序列和原始ID
    用于alignment_quality策略
    性能优化：批量处理序列
    
    Args:
        file_path: FASTA文件路径
        blast_params: BLAST参数
    
    Returns:
        物种到ID-序列字典的映射
    """
    species_to_seqs = defaultdict(dict)
    
    try:
        # 性能优化：使用生成器并批量处理
        record_count = 0
        batch_size = 100
        batch = []
        
        for record in SeqIO.parse(file_path, "fasta"):
            batch.append(record)
            record_count += 1
            
            if len(batch) >= batch_size:
                process_duplicate_batch(batch, species_to_seqs, file_path, blast_params)
                batch = []
        
        # 处理最后一批
        if batch:
            process_duplicate_batch(batch, species_to_seqs, file_path, blast_params)
            
        logger.debug(f"在 {file_path} 中找到并处理了 {record_count} 个序列")
        
        # 记录发现的重复
        duplicates = {species: len(ids) for species, ids in species_to_seqs.items() if len(ids) > 1}
        if duplicates:
            logger.info(f"在 {file_path} 中发现重复物种: {duplicates}")
            
        return species_to_seqs
    except Exception as e:
        logger.error(f"解析 FASTA 文件 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return {}

def process_duplicate_batch(batch, species_to_seqs, file_path, blast_params):
    """
    处理包含可能重复的FASTA记录批次
    
    Args:
        batch: 要处理的记录批次
        species_to_seqs: 物种到ID-序列字典的映射
        file_path: 原始文件路径
        blast_params: BLAST参数
    """
    for record in batch:
        # 从头部提取物种名
        species = extract_species_name(record.description)
        record_id = record.id
        sequence = str(record.seq)
        
        # 预处理序列
        sequence, _ = preprocess_cds_sequence(sequence, record_id, file_path, blast_params)
        
        # 使用原始ID存储以跟踪它们
        species_to_seqs[species][record_id] = sequence

def safe_translate_codon(codon):
    """
    安全地将密码子翻译为氨基酸，处理模糊碱基
    
    Args:
        codon: 三字母核苷酸密码子
        
    Returns:
        单字母氨基酸或'X'表示未知
    """
    # 转为大写并移除空白
    codon = codon.upper().strip()
    
    # 检查密码子是否存在于查找表中
    if codon in GENETIC_CODE:
        return GENETIC_CODE[codon]
    
    # 如果密码子包含N或其他模糊核苷酸
    if 'N' in codon or any(base not in 'ACGT' for base in codon):
        return 'X'
        
    # 这种情况不应该出现在有效的DNA序列中
    logger.warning(f"未知的密码子: '{codon}'")
    return 'X'

def convert_duplicates_to_single_sequences(species_to_seqs, strategy='longest'):
    """
    将包含重复序列的字典转换为每个物种一个序列的字典
    
    Args:
        species_to_seqs: 物种到ID-序列字典的映射
        strategy: 选择策略 ('longest', 'first')
    
    Returns:
        物种到单个序列的字典
    """
    single_seqs = {}
    
    for species, id_to_seq in species_to_seqs.items():
        if not id_to_seq:
            continue
            
        if len(id_to_seq) == 1:
            # 只有一个序列
            single_seqs[species] = next(iter(id_to_seq.values()))
        else:
            # 多个序列，根据策略选择
            if strategy == 'longest':
                # 选择最长的序列
                best_id = max(id_to_seq.items(), key=lambda x: len(x[1]))[0]
                single_seqs[species] = id_to_seq[best_id]
                logger.debug(f"为物种 {species} 选择最长序列 {best_id} (长度: {len(id_to_seq[best_id])})")
            elif strategy == 'first':
                # 选择第一个序列
                first_id = next(iter(id_to_seq))
                single_seqs[species] = id_to_seq[first_id]
                logger.debug(f"为物种 {species} 选择第一个序列 {first_id}")
            else:
                # 默认使用最长策略
                best_id = max(id_to_seq.items(), key=lambda x: len(x[1]))[0]
                single_seqs[species] = id_to_seq[best_id]
                logger.debug(f"为物种 {species} 使用默认策略选择序列 {best_id}")
    
    return single_seqs

def write_output(species_seqs, output_file):
    """
    将处理后的序列写入FASTA文件
    修复了处理字典类型的问题
    性能优化：批量写入
    """
    if not species_seqs:
        logger.warning(f"没有序列可写入 {output_file}")
        return False
        
    try:
        # 如果不存在则创建目录
        os.makedirs(os.path.dirname(os.path.abspath(output_file)), exist_ok=True)
        
        # 检查species_seqs的格式，如果是嵌套字典则转换
        if species_seqs and isinstance(next(iter(species_seqs.values())), dict):
            # 这是alignment_quality策略返回的格式，需要转换
            logger.debug("检测到嵌套字典格式，转换为单序列格式")
            species_seqs = convert_duplicates_to_single_sequences(species_seqs, strategy='longest')
        
        # 批量写入以提高性能
        batch_size = 10  # 每批处理的物种数
        species_list = sorted(species_seqs.keys())
        
        with open(output_file, 'w') as f:
            for i in range(0, len(species_list), batch_size):
                batch = species_list[i:i+batch_size]
                
                for species in batch:
                    seq = species_seqs[species]
                    # 确保seq是字符串
                    if not isinstance(seq, str):
                        logger.error(f"物种 {species} 的序列不是字符串类型: {type(seq)}")
                        continue
                        
                    f.write(f">{species}\n")
                    # 以60个字符一行写入序列，提高可读性
                    for j in range(0, len(seq), 60):
                        f.write(seq[j:j+60] + '\n')
                
                # 定期刷新缓冲区
                f.flush()
        
        seq_length = len(next(iter(species_seqs.values()))) if species_seqs else 0
        logger.info(f"成功写入 {len(species_seqs)} 条序列到 {output_file} (长度: {seq_length} bp)")
        return True
    except Exception as e:
        logger.error(f"写入 {output_file} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return False

def process_cds_files(input_dir, output_dir, duplicate_strategy='longest', 
                     blast_params=None, threads=4):
    """
    处理目录中的所有CDS文件，执行序列预处理和重复物种处理
    
    Args:
        input_dir: 包含CDS文件的输入目录
        output_dir: 输出目录
        duplicate_strategy: 重复物种处理策略
        blast_params: BLAST参数字典
        threads: 并行处理线程数
        
    Returns:
        处理结果字典
    """
    # 查找所有CDS文件
    cds_files = find_cds_files(input_dir)
    if not cds_files:
        logger.error("未找到.cds文件。退出。")
        return {}
        
    logger.info(f"将处理 {len(cds_files)} 个CDS文件")
    
    # 创建输出目录
    processed_dir = os.path.join(output_dir, "processed")
    os.makedirs(processed_dir, exist_ok=True)
    
    results = {}
    
    # 使用线程池并行处理每个CDS文件
    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = {}
        
        for file_name in cds_files:
            file_path = os.path.join(input_dir, file_name)
            gene_name = file_name.replace('.cds', '')
            output_file = os.path.join(processed_dir, f"{gene_name}_processed.fasta")
            
            futures[executor.submit(
                process_single_cds_file,
                file_path,
                output_file,
                duplicate_strategy,
                blast_params
            )] = file_name
        
        # 收集结果并显示进度
        completed = 0
        total_futures = len(futures)
        
        for future in futures:
            file_name = futures[future]
            try:
                result = future.result()
                completed += 1
                progress = completed / total_futures * 100
                logger.info(f"进度: {completed}/{total_futures} 文件处理 ({progress:.1f}%)")
                
                if result:
                    results[file_name] = result
                    logger.info(f"处理 {file_name}: {result['species_count']} 个物种")
                else:
                    logger.error(f"处理 {file_name} 失败")
            except Exception as e:
                logger.exception(f"处理 {file_name} 时出错: {str(e)}")
    
    return results

def process_single_cds_file(file_path, output_file, duplicate_strategy, blast_params):
    """
    处理单个CDS文件
    
    Args:
        file_path: CDS文件路径
        output_file: 输出文件路径
        duplicate_strategy: 重复处理策略
        blast_params: BLAST参数
        
    Returns:
        处理结果字典或None
    """
    try:
        file_name = os.path.basename(file_path)
        logger.debug(f"开始处理 {file_name}")
        
        # 解析FASTA文件
        species_seqs = parse_fasta(file_path, duplicate_strategy, blast_params)
        
        if not species_seqs:
            logger.error(f"无法从 {file_path} 解析序列")
            return None
            
        logger.info(f"从 {file_name} 解析得到 {len(species_seqs)} 个物种的序列")
        
        # 写入处理后的序列
        success = write_output(species_seqs, output_file)
        if not success:
            logger.error(f"写入输出文件 {output_file} 失败")
            return None
        
        # 返回结果
        result = {
            'input_file': file_path,
            'output_file': output_file,
            'species_count': len(species_seqs),
            'sequences': species_seqs
        }
        
        return result
        
    except Exception as e:
        logger.error(f"处理 {file_path} 时出错: {str(e)}")
        logger.error(traceback.format_exc())
        return None

def main():
    """主函数：解析参数并执行序列预处理"""
    parser = argparse.ArgumentParser(description="CDS序列预处理和重复物种处理")
    parser.add_argument("--input_dir", default=".", help="包含CDS文件的目录")
    parser.add_argument("--output_dir", default="./output", help="输出文件目录")
    parser.add_argument("--threads", type=int, default=4, help="用于并行处理的线程数")
    
    # BLAST参数
    parser.add_argument("--use_blast", action="store_true", help="使用BLAST帮助修复非3倍数序列")
    parser.add_argument("--blast_path", help="BLASTN可执行文件的绝对路径")
    parser.add_argument("--blast_db_cache_size", type=int, default=5, help="BLAST数据库缓存大小")
    
    # 重复处理策略
    parser.add_argument("--duplicate_strategy", 
                      choices=['longest', 'first', 'rename', 'alignment_quality'], 
                      default='longest',
                      help="处理文件中重复物种的策略")
    
    # 性能优化参数
    parser.add_argument("--batch_size", type=int, default=100, help="批处理大小，用于优化内存使用")
    parser.add_argument("--memory_limit", type=int, default=0, 
                      help="内存限制（MB），0表示无限制。当处理大文件时有用。")
    
    # 其他参数
    parser.add_argument("--log_level", choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO',
                      help="设置日志级别")
    parser.add_argument("--clean_temp", action="store_true", default=True,
                      help="处理后清理临时文件")
    parser.add_argument("--min_repeat_report", type=int, default=70,
                      help="报告为质量问题所需的最小重复数量 (默认: 70)")
    
    args = parser.parse_args()
    
    # 更新基于命令行参数的全局常量
    global MIN_REPEAT_COUNT_TO_REPORT, BLAST_DB_CACHE_SIZE
    MIN_REPEAT_COUNT_TO_REPORT = args.min_repeat_report
    BLAST_DB_CACHE_SIZE = args.blast_db_cache_size
    
    # 如果启用，验证BLAST
    if args.use_blast and not args.blast_path:
        parser.error("启用BLAST处理时，必须提供--blast_path参数")
    
    # 设置日志级别
    logging.getLogger().setLevel(getattr(logging, args.log_level))
    
    # 配置内存限制（如果指定）
    if args.memory_limit > 0:
        import resource
        # 将软限制设置为用户请求的限制（以字节为单位）
        resource.setrlimit(
            resource.RLIMIT_AS, 
            (args.memory_limit * 1024 * 1024, resource.RLIM_INFINITY)
        )
        logger.info(f"设置内存限制为 {args.memory_limit} MB")
    
    # 记录启动时间和参数
    logger.info(f"启动CDS序列预处理，参数: {vars(args)}")
    start_time = time.time()
   
    # 设置目录结构
    try:
        # 确保输出目录存在
        os.makedirs(args.output_dir, exist_ok=True)
        
        # 创建必要的子目录
        subdirs = ["processed", "temp"]
        for subdir in subdirs:
            full_path = os.path.join(args.output_dir, subdir)
            os.makedirs(full_path, exist_ok=True)
            logger.debug(f"创建目录: {full_path}")
    except Exception as e:
        logger.error(f"设置输出目录时出错: {str(e)}")
        return 1
    
    # 准备BLAST参数
    blast_params = None
    if args.use_blast:
        blast_params = {
            'use_blast': args.use_blast,
            'blast_path': args.blast_path,
            'temp_dir': os.path.join(args.output_dir, "temp")
        }
    
    # 处理所有CDS文件
    results = process_cds_files(
        args.input_dir,
        args.output_dir,
        args.duplicate_strategy,
        blast_params,
        args.threads
    )
    
    if not results:
        logger.error("没有文件成功处理")
        return 1
    
    # 如果请求，清理临时文件
    if args.clean_temp:
        try:
            logger.info("清理临时文件")
            temp_dir = os.path.join(args.output_dir, "temp")
            if os.path.exists(temp_dir):
                shutil.rmtree(temp_dir)
        except Exception as e:
            logger.warning(f"清理临时文件时出错: {str(e)}")
    
    # 记录执行时间
    end_time = time.time()
    execution_time = end_time - start_time
    hours, remainder = divmod(execution_time, 3600)
    minutes, seconds = divmod(remainder, 60)
    
    logger.info(f"预处理完成。耗时: {int(hours)}小时 {int(minutes)}分钟 {seconds:.1f}秒")
    logger.info(f"成功处理 {len(results)} 个文件")
    logger.info(f"结果保存在: {args.output_dir}")
    
    return 0

if __name__ == "__main__":
    try:
        # 启用垃圾收集诊断，帮助调试内存问题
        import gc
        gc.set_debug(gc.DEBUG_STATS)
        
        sys.exit(main())
    except Exception as e:
        logger.critical(f"未捕获的异常: {str(e)}")
        logger.critical(traceback.format_exc())
        sys.exit(1) 
