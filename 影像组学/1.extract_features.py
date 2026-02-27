"""
影像组学特征提取脚本

从 NIfTI 格式的 CT 影像和分割掩膜中，使用 PyRadiomics 批量提取影像组学特征，
支持多进程并行加速，结果保存为 CSV 文件。

用法:
    python extract_features.py

    通过环境变量 RADIOMICS_N_JOBS 控制并行进程数（默认 6，设为 1 禁用并行）:
        RADIOMICS_N_JOBS=1 python extract_features.py
"""

import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import pandas as pd
import yaml
from radiomics import featureextractor

# ========== 1. 全局配置

DATA_ROOT = Path("data/nifti")  # NIfTI 数据根目录
IMAGE_FILENAME = "ct.nii.gz"  # CT文件名
MASK_FILENAME = "seg.nii.gz"  # 分割文件名
PARAMS_PATH = Path("params.yaml")  # 参数文件
OUTPUT_PATH = Path("输出/features.csv")  # 输出CSV

# 并行进程数，设为 1 时禁用，根据 CPU 核心数和内存调整，设置过大小心卡死系统
N_JOBS = int(os.getenv("RADIOMICS_N_JOBS", "6"))

# ========== 2. 类型别名

# 单个提取任务: (case_id, image_path, mask_path)
Task = tuple[str, str, str]

# image 与 mask 配对，注意，必须提前清洗好数据，保证每个子文件夹内都有配对的 ct.nii.gz 和 seg.nii.gz 两个文件
ImageMaskPair = tuple[Path, Path]

# ========== 3. 全局状态

# 全局特征提取器实例，子进程初始化时创建，供 extract_one_case 使用
_EXTRACTOR = None

# ========== 4. 数据发现与 ID 生成

def collect_image_mask_pairs(data_root: Path) -> list[ImageMaskPair]:
    """递归扫描 data_root，收集所有 (影像, 掩膜) 文件对。"""

    # pairs本身为列表，元素为CT与分割文件的配对
    pairs: list[ImageMaskPair] = []

    for image_path in sorted(data_root.rglob(IMAGE_FILENAME)):

        mask_path = image_path.with_name(MASK_FILENAME)

        # 配对CT与分割文件
        pairs.append((image_path, mask_path))

    return pairs


def make_case_id(data_root: Path, image_path: Path) -> str:
    """根据影像路径相对于 data_root 的位置生成病例 ID。

    例: data_root/patient01/study01/ct.nii.gz → "patient01/study01"
    """
    return image_path.parent.relative_to(data_root).as_posix()

# ========== 5. 特征提取

def _init_worker(params_path: str) -> None:
    """子进程初始化函数：加载参数并创建特征提取器实例。

    每个子进程只调用一次，避免重复初始化的开销。
    """

    global _EXTRACTOR

    # 读取yaml参数文件，创建特征提取器实例
    with Path(params_path).open("r", encoding="utf-8") as f:
        params = yaml.safe_load(f)

    _EXTRACTOR = featureextractor.RadiomicsFeatureExtractor(params)


def extract_one_case(task: Task) -> dict:
    """对单个病例执行特征提取。

    Args:
        task: (case_id, image_path, mask_path) 三元组。

    Returns:
        成功时返回 {"ok": True, "row": {特征字典}}，
        失败时返回 {"ok": False, "case_id": ..., "error": ...}。
    """

    case_id, image_path, mask_path = task

    # 确保提取器已初始化（子进程中调用 _init_worker 后才会有值）
    # 防御性编程，理论上不应该发生，但是Pylance会警告未赋值的变量，最好加上
    assert _EXTRACTOR is not None
    try:
        result = _EXTRACTOR.execute(image_path, mask_path)
    except Exception as exc:
        return {"ok": False, "case_id": case_id, "error": str(exc)}

    # 构建结果行，包含病例ID、文件路径和提取的特征
    row = {"case_id": case_id, "image_path": image_path, "mask_path": mask_path}
    row.update(result)
    
    # 成功返回结果
    return {"ok": True, "row": row}

# ========== 6. 结果收集与日志

def _process_output(
    output: dict, records: list[dict], idx: int, total: int
) -> None:
    """处理单条提取结果：成功则收集，失败则打印错误。"""
    if output["ok"]:
        records.append(output["row"])
    else:
        print(f"  [{idx}/{total}] [error] {output['case_id']}: {output['error']}")

# ========== 7. 串行 / 并行执行

def _run_serial(tasks: list[Task]) -> list[dict]:
    """串行单进程：逐个提取特征（N_JOBS=1 时使用）。"""

    _init_worker(str(PARAMS_PATH))
    records: list[dict] = []

    for idx, task in enumerate(tasks, start=1):
        output = extract_one_case(task)
        _process_output(output, records, idx, len(tasks))
    
    return records


def _run_parallel(tasks: list[Task]) -> list[dict]:
    """并行加速：使用多进程池加速特征提取。"""

    records: list[dict] = []

    with ProcessPoolExecutor(max_workers=N_JOBS, initializer=_init_worker, initargs=(str(PARAMS_PATH),),) as executor:

        for idx, output in enumerate(executor.map(extract_one_case, tasks), start=1):
            case_id = tasks[idx - 1][0]
            _process_output(output, records, idx, len(tasks))

    return records

# ========== 8. 结果保存

def _save_features(records: list[dict]) -> None:
    """将提取结果整理为 DataFrame 并保存为 CSV。

    列顺序: case_id, image_path, mask_path, ...特征列...
    """

    df = pd.DataFrame(records)
    metadata_cols = ["case_id", "image_path", "mask_path"]

    # 将特征列放在 metadata 列之后，保持 CSV 可读性
    feature_cols = [col for col in df.columns if col not in metadata_cols]
    df = df[metadata_cols + feature_cols]

    df.to_csv(OUTPUT_PATH, index=False)

    print(f"共 {len(df)} 条记录已保存至: {OUTPUT_PATH}")


# ========== 0. 主入口


def main() -> None:
    """主入口：收集数据 → 提取特征 → 保存结果。"""

    # 1. 收集影像/掩膜配对
    pairs = collect_image_mask_pairs(DATA_ROOT)

    tasks: list[Task] = [
        (make_case_id(DATA_ROOT, img), str(img), str(msk))
        for img, msk in pairs
    ]

    # 2. 执行特征提取，根据 N_JOBS 决定串行单进程还是并行加速
    print(f"对 {len(tasks)} 对影像/掩膜对，使用 {N_JOBS} 个进程")
    records = _run_serial(tasks) if N_JOBS <= 1 else _run_parallel(tasks)

    # 3. 保存结果
    _save_features(records)

if __name__ == "__main__":
    main()
