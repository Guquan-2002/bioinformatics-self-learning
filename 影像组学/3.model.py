import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split, StratifiedKFold, cross_val_score
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import warnings
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score,
    f1_score, roc_auc_score, classification_report,
    roc_curve, auc, confusion_matrix
)
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.family'] = 'LXGW WenKai Mono'
matplotlib.rcParams['font.size'] = 12

warnings.filterwarnings('ignore')

# ========== 1. 加载数据 ==========
# 真实数据直接导入CSV文件，需要包含标签和特征
# df = pd.read_csv("features_selected.csv")

# 模拟数据演示（1000个样本，10个特征，二分类）
np.random.seed(412)
n_samples = 1000
n_features = 10

# 模拟特征数据
X = pd.DataFrame(
    np.random.randn(n_samples, n_features),
    columns=[f"feature_{i}" for i in range(n_features)]
)

# 模拟二分类标签
y = pd.Series(np.random.randint(0, 2, size=n_samples), name="label")

print(f"样本数: {len(X)}, 特征数: {X.shape[1]}")
print(f"标签分布:\n{y.value_counts()}")

# ========== 2. 数据集划分 ==========

X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.3, random_state=412, stratify=y)

print(f"\n训练集: {len(X_train)}, 测试集: {len(X_test)}")

# ========== 3. 定义模型 ==========
"""
Logistic Regression 可解释性强
SVM 高维小样本表现好
Random Forest 能输出特征重要性
XGBoost 小样本容易过拟合，大数据集表现好
LightGBM  比 XGBoost 快，大数据集优势明显
KNN 简单直观，但对特征尺度敏感
Naive Bayes 需要假设特征独立
Decision Tree 单棵树容易过拟合
AdaBoost 对噪声敏感
MLP 需要较多样本
"""
models = {
    # 逻辑回归，可解释性强，但只能处理线性关系
    "Logistic Regression": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", LogisticRegression(max_iter=1000, random_state=412))
    ]),
    
    # 随机森林，能处理非线性关系
    "Random Forest": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", RandomForestClassifier(n_estimators=100, random_state=412))
    ]),

    # 支持向量机，泛化能力强，可解释性较差
    "SVM": Pipeline([
        ("scaler", StandardScaler()),
        ("clf", SVC(kernel="rbf", probability=True, random_state=412))
    ]),
}

# ========== 4. 多个模型交叉验证 ==========
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=412)
cv_results = {}

for name, model in models.items():
    scores = cross_val_score(model, X_train, y_train, cv=cv, scoring="roc_auc")
    cv_results[name] = scores
    print(f"{name}:")
    print(f"  AUC = {scores.mean():.3f} ± {scores.std():.3f}")

# ========== 5. 训练最终模型 ==========
for name, model in models.items():
    model.fit(X_train, y_train)

best_model_name = max(cv_results, key=lambda k: cv_results[k].mean())
best_model = models[best_model_name]

print(f"\n最佳模型: {best_model_name}")

# =============================================
# 图1：ROC 曲线
# =============================================
colors = {"Logistic Regression": "#66ccff", "Random Forest": "#FFCC66", "SVM": "#ee0000"}

fig, axes = plt.subplots(1, 2, figsize=(14, 6))

for ax, (X_eval, y_eval, title) in zip(axes, [
    (X_train, y_train, "训练集"),
    (X_test, y_test, "测试集"),
]):
    for name, model in models.items():
        y_prob_eval = model.predict_proba(X_eval)[:, 1]
        fpr, tpr, _ = roc_curve(y_eval, y_prob_eval)
        roc_auc = auc(fpr, tpr)
        ax.plot(fpr, tpr, linewidth=2, color=colors[name],
                label=f"{name} (AUC={roc_auc:.3f})")

    ax.plot([0, 1], [0, 1], 'k--', linewidth=1, alpha=0.5)
    ax.set_xlabel("1 - 特异度 (FPR)")
    ax.set_ylabel("灵敏度 (TPR)")
    ax.set_title(f"{title} ROC 曲线")
    ax.legend(loc="lower right", fontsize=10)
    ax.set_xlim([-0.02, 1.02])
    ax.set_ylim([-0.02, 1.02])

plt.tight_layout()
plt.show()

# =============================================
# 图2：混淆矩阵
# =============================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

for idx, (name, model) in enumerate(models.items()):
    y_pred_m = model.predict(X_test)
    cm = confusion_matrix(y_test, y_pred_m)
    ax = axes[idx]
    im = ax.imshow(cm, interpolation='nearest', cmap="Blues")

    for i in range(2):
        for j in range(2):
            color = "#EE0000" if cm[i, j] > cm.max() / 2 else "black"
            ax.text(j, i, str(cm[i, j]), ha="center", va="center",
                    fontsize=20, fontweight="bold", color=color)

    ax.set_title(name, fontsize=13)
    ax.set_xlabel("预测值")
    ax.set_ylabel("真实值")
    ax.set_xticks([0, 1])
    ax.set_yticks([0, 1])
    ax.set_xticklabels(["阴性", "阳性"])
    ax.set_yticklabels(["阴性", "阳性"])

plt.tight_layout()
plt.show()

# =============================================
# 图3：多模型性能对比
# =============================================
metrics_data = []
for name, model in models.items():
    y_pred_m = model.predict(X_test)
    y_prob_m = model.predict_proba(X_test)[:, 1]
    metrics_data.append({
        "Model": name,
        "正确率": accuracy_score(y_test, y_pred_m),
        "精确率": precision_score(y_test, y_pred_m),
        "召回率": recall_score(y_test, y_pred_m),
        "F1": f1_score(y_test, y_pred_m),
        "AUC": roc_auc_score(y_test, y_prob_m),
    })

df_metrics = pd.DataFrame(metrics_data)
print("\n===== 各模型测试集指标 =====")
print(df_metrics.to_string(index=False))

metric_names = ["正确率", "精确率", "召回率", "F1", "AUC"]
x = np.arange(len(metric_names))
width = 0.25

fig, ax = plt.subplots(figsize=(12, 6))
for i, (_, row) in enumerate(df_metrics.iterrows()):
    values = [row[m] for m in metric_names]
    bars = ax.bar(x + i * width, values, width,
                  label=row["Model"], color=list(colors.values())[i])
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                f"{val:.2f}", ha='center', va='bottom', fontsize=9)

ax.set_ylabel("得分")
ax.set_title("各模型性能对比")
ax.set_xticks(x + width)
ax.set_xticklabels(metric_names)
ax.set_ylim(0, 1.15)
ax.legend()
ax.grid(axis='y', alpha=0.3)

plt.tight_layout()
plt.show()

# =============================================
# 图4：随机森林特征重要性
# =============================================
rf_clf = models["Random Forest"].named_steps["clf"]
importances = rf_clf.feature_importances_
indices = np.argsort(importances)[::-1]

fig, ax = plt.subplots(figsize=(10, 6))
bars = ax.barh(range(len(indices)), importances[indices], color="#66CCFF")
ax.set_yticks(range(len(indices)))
ax.set_yticklabels([X.columns[i] for i in indices])
ax.set_xlabel("重要性")
ax.set_title("随机森林特征重要性")
ax.invert_yaxis()

# 标数值
for bar, val in zip(bars, importances[indices]):
    ax.text(bar.get_width() + 0.002, bar.get_y() + bar.get_height() / 2,
            f"{val:.3f}", va='center', fontsize=10)

plt.tight_layout()
plt.show()