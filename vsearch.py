import os
import subprocess
from glob import glob
from concurrent.futures import ThreadPoolExecutor
import multiprocessing

# ======================
# 硬编码配置参数
# ======================
INPUT_DIR = "locus_fasta_files"  # 输入目录
OUTPUT_DIR = "locus_cluster_results"  # 输出目录
VSEARCH_PATH = "/share/home/zhanglab/user/yangting/software/vsearch/bin/vsearch"  # VSEARCH可执行文件路径
MAX_WORKERS = 4  # 并行线程数
IDENTITY_THRESHOLD = 0.9  # 序列相似性阈值


def run_vsearch_on_locus(fasta_file, output_dir, threads=1, identity=0.85):
    """对单个位点的FASTA文件执行VSEARCH聚类"""
    VSEARCH_PATH = r"F:\华大基因\tools\vsearch-2.30.0-win-x86_64\vsearch-2.30.0-win-x86_64\bin\vsearch.exe"

    locus_id = os.path.basename(fasta_file).replace('.fa', '')
    output_prefix = os.path.join(output_dir, locus_id)

    # 构建VSEARCH命令
    command = [
        VSEARCH_PATH,  # 指定完整路径+可执行文件名
        "--cluster_fast", fasta_file,
        "--id", str(identity),
        "--strand", "both",
        "--uc", f"{output_prefix}.uc",
        "--threads", str(threads),
        "--maxseqlength", "100000",
        "--minseqlength", "2",
        "--sizeout",
    ]

    print(f"Running command: {' '.join(command)}")

    try:
        # 使用subprocess代替os.system，便于错误处理
        result = subprocess.run(
            command,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print(f"✓ 完成位点: {locus_id}")
        return True
    except subprocess.CalledProcessError as e:
        print(f"✗ 位点 {locus_id} 聚类失败: {e.stderr}")
        return False


def batch_run_vsearch(input_dir, output_dir, max_workers, identity=0.85):
    """批量处理所有位点"""
    os.makedirs(output_dir, exist_ok=True)
    fasta_files = glob(os.path.join(input_dir, "*.fa"))
    total_files = len(fasta_files)

    if total_files == 0:
        print(f"错误: 输入目录 {input_dir} 中无FASTA文件")
        return

    print(f"开始处理 {total_files} 个位点...")
    threads_per_task = max(1, multiprocessing.cpu_count() // max_workers)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(run_vsearch_on_locus, fasta_file, output_dir, threads_per_task, identity): fasta_file
            for fasta_file in fasta_files
        }

        completed_count = 0
        for future in futures:
            future.result()
            completed_count += 1
            print(f"进度: {completed_count}/{total_files}")

    print(f"\n所有位点处理完成！结果保存至: {output_dir}")


if __name__ == "__main__":
    # 运行批量聚类
    batch_run_vsearch(
        input_dir=INPUT_DIR,
        output_dir=OUTPUT_DIR,
        max_workers=MAX_WORKERS,
        identity=IDENTITY_THRESHOLD
    )