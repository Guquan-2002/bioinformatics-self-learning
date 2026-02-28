import pandas as pd
import numpy as np
from scipy import stats
from sklearn.preprocessing import StandardScaler

def split_cols(df: pd.DataFrame):
    """将列分为元数据列、诊断列和特征列，返回特征子集。"""

    meta_cols = ["case_id", "image_path", "mask_path"]
    diag_cols = [col for col in df.columns if col.startswith("diagnostics_")]
    feature_cols = [col for col in df.columns if col not in meta_cols + diag_cols]

    print(f"原始列数: {len(df.columns)}")
    print(f"诊断信息列: {len(diag_cols)}")
    print(f"特征列数: {len(feature_cols)}")
    print(f"样本数: {len(df)}")
    
    return meta_cols, diag_cols, feature_cols

def missing(df: pd.DataFrame):
    """检查是否存在缺失值的特征数。"""

    missing_cols = df.columns[df.isnull().any()].tolist()

    print(f"\n有缺失值的特征数: {len(missing_cols)}")
    if missing_cols:
        print(missing_cols)

    return missing_cols

def constant(df: pd.DataFrame):
    """检查是否存在常量特征。"""

    constant_cols = df.columns[df.nunique() <= 1].tolist()

    print(f"\n常量特征数: {len(constant_cols)}")
    if constant_cols:
        print(constant_cols)
    
    return constant_cols

def infinite(df: pd.DataFrame):
    """检查是否存在无穷值。"""

    inf_mask = np.isinf(df.select_dtypes(include=[np.number])).any()
    inf_cols = inf_mask[inf_mask].index.tolist()

    print(f"\n含无穷值的特征数: {len(inf_cols)}")
    if inf_cols:
        print(inf_cols)

    return inf_cols

def replace_outliers(df: pd.DataFrame) -> pd.DataFrame:
    """将 Z-score 绝对值超过 3 的异常值替换为该列中位数。"""

    replaced_count = 0
    for col in df.columns:
        col_z = np.abs(stats.zscore(df[col], nan_policy='omit'))
        median_val = df[col].median()
        mask = col_z > 3
        replaced_count += mask.sum()
        df.loc[mask, col] = median_val
    print(f"\n|Z| ＞ 3的特征数: {replaced_count} ")

    return df

def standardize(df: pd.DataFrame) -> pd.DataFrame:
    """Z-score 标准化，使每个特征均值≈0，标准差≈1。"""

    scaler = StandardScaler()
    return pd.DataFrame(
        scaler.fit_transform(df),
        columns=df.columns,
        index=df.index
    )

def main():
    """主入口"""

    # 读取数据
    df = pd.read_csv("输出/features.csv")

    # 分离元数据和特征，只保留特征列进行检查
    meta_cols, diag_cols, feature_cols = split_cols(df)

    # copy()创建特征子集进行检查，避免修改原始 DataFrame
    df_features: pd.DataFrame = df[feature_cols].copy()  # type: ignore[assignment]

    # 删除不良特征，此处包括缺失值、常量、无穷值特征
    df_features = df_features.drop(columns=missing(df_features))
    df_features = df_features.drop(columns=constant(df_features))
    df_features = df_features.drop(columns=infinite(df_features))

    # 替换异常值，使用 Z-score 方法，|Z| ＞ 3 的视为异常值，使用中位值替换
    df_features = replace_outliers(df_features)

    # 标准化后每个特征均值≈0，标准差≈1
    df_scaled = standardize(df_features)

    # 把 case_id 加回来，方便后续追溯
    df_clean = pd.concat([df[["case_id"]], df_scaled], axis=1)
    df_clean.to_csv("输出/features_cleaned.csv", index=False)
    
    print(f"\n最终: {len(df_clean)} 个样本, {len(df_scaled.columns)} 个特征")

if __name__ == "__main__":
    main()