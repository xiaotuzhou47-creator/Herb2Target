# Herb2Target

面向中药复方网络药理学的端到端集成计算平台。

Herb2Target 将草药目录、ADME 筛选、靶点预测、PPI 网络、口袋对接和 3D 可视化整合到一个工作流中，支持单味药、经典方剂和自定义 SMILES 输入。

## 系统要求

- Python 3.10+
- Flask（Web 框架）
- RDKit（化学信息学）
- AutoDock Vina 1.1.2（分子对接）
- OpenBabel 3.x（配体准备）

## 快速开始

### 安装依赖

```bash
pip install flask flask-cors rdkit openbabel
```

AutoDock Vina 1.1.2 需另行安装，详见 http://vina.scripps.edu/

### 运行

```bash
python app.py
```

启动后访问 http://localhost:5000

## 项目结构

| 文件 | 说明 |
|------|------|
| `app.py` | Web 应用主程序（Flask） |
| `docking.py` | 分子对接引擎（AutoDock Vina） |
| `tcm_demo.db` | SQLite 数据库（约 2.2 MB） |
| `itmap_raw.json` | 参考化合物-靶点映射（765 化合物，72 靶点） |
| `itmap_stratified.json` | 按化学类型分层的参考库 |
| `itmap.json` | 简化的化合物-靶点映射 |
| `run_validation.py` | 10 个验证化合物的靶点预测重现脚本 |
| `templates/index.html` | Web 前端页面 |
| `static/` | 前端资源（Vue 3, Element Plus, axios） |

## 数据库内容

- 草药: 113 味（TCMSP + HERB 来源）
- 活性化合物: 920 个（经 ADME 筛选: OB >= 30%, DL >= 0.18）
- 方剂: 121 首（15 个功效类别）
- 参考库: 765 个化合物，72 个药物靶点（ChEMBL + NPASS 来源）
- 民族药: 14 种苗族药材

## 靶点预测算法

1. 查询化合物编码为 2048 位 ECFP4 指纹
2. 与 765 个参考化合物计算 Tanimoto 相似度
3. Top-30 相似参考化合物按相似度平方加权投票
4. TF-IDF 风格归一化（除 log₂(频次)）
5. Fpocket 可药性评分（PS）加权的综合排序
6. 返回 Top-15 靶点

## 验证

```bash
python run_validation.py
```

重现论文中 10 个天然产物-靶点对的验证结果。

## 补充材料

论文补充表 S1-S7 提供了以下详细信息：
- S1: 377 化合物 × 13 靶点对接基准数据集
- S2: 753 化合物参考库
- S3: 跨草药验证
- S4: 验证化合物-靶点对详情
- S5: 结构生物学合理性论证
- S6: SwissTargetPrediction 对比
- S7: 完整方剂目录

## 注意事项

本应用仅为演示与验证用途，部分模板文件（FinanceGen、PptxGen、ProjectPlanGen）与核心功能无关。

## 许可证

MIT License
