# HapBlock V17 中文使用手册

版本：17.0.0  
日期：2026-10-20

## 1. 流程用途

HapBlock V17 从群体 VCF 中检测栽培群体的基因组 haplotype block，并回答以下问题：

1. 每个 block 位于哪里、长度是多少；
2. 每个 block 中有多少个主要 haplotype；
3. 每个 haplotype 的频率及其携带样本；
4. 几 Mb 到几十 Mb 的大 block 内是否还有更小的子 block；
5. 哪些 block 满足后续性状关联分析的样本比例要求；
6. 提供 wild 样本时，哪些 wild 样本与哪个 cultivar haplotype 最接近。

流程支持两种自动模式：

- `CULTIVAR_ONLY`：只给 cultivar list，检测 block、haplotype、近期亚型和关联分析矩阵。
- `CULTIVAR_WILD`：同时给 cultivar list 与 wild list。block 和 cultivar haplotype 仍只由 cultivar 决定，完成后才投影 wild。

这里的 `Hap01` 等名称表示一个区段内相似的二倍体基因型模式组。对于高纯合的自交材料，这通常可以近似表示单倍型组。V17 当前不拆分未相位杂合基因型中的两条染色体，因此不应把结果解释为实验相位得到的单条染色体序列。

## 2. V17 的主要改变

相对于旧流程，V17 的关键改变是：

- 使用 complete-linkage（完全连接）聚类，保证同组最远样本对也要通过距离阈值；
- 避免 `A` 接近 `B`、`B` 接近 `C`，但 `A` 不接近 `C` 时仍把三者连成一组；
- 不再限制最多只能报告几个 haplotype；
- 一次读取 VCF 并建立紧凑缓存，后续窗口不再反复解压整个 VCF；
- 同时保留 coarse 父 block 和 fine 子 block，fine 扫描最低达到 10 kb；
- 短的低信息区可作为内部 gap 保留，不会轻易切断几十 Mb 的大 block；
- common SNP 定义较古老的 backbone haplotype，再在每个 backbone 内使用组内多态 SNP 检测 recent subhaplotype；
- 直接输出可用于性状关联的因子矩阵和 one-hot 矩阵。

## 3. 安装

推荐使用 Conda：

```bash
tar -xzf hapblock_v17.0.0.tar.gz
cd hapblock_v17_release
conda env create -f environment.yml
conda activate hapblock-v17
python3 hapblock_v17.py --version
```

也可以使用 Python virtual environment：

```bash
python3 -m venv .venv
source .venv/bin/activate
python3 -m pip install -r requirements.txt
python3 hapblock_v17.py --version
```

要求 Python 3.10 或更高版本。核心依赖为 NumPy、SciPy、scikit-learn 和 Matplotlib。

## 4. 输入文件

### 4.1 VCF

输入应为按染色体和坐标排序的普通 VCF 或 `.vcf.gz` 文件。建议：

- 使用未经过 MAF 过滤的群体 VCF；
- 至少包含 `GT` 字段；
- cultivar 和可选 wild 样本在同一个 VCF 中；
- VCF header 中包含 `##contig=<ID=...,length=...>`；
- V17 当前保留双等位 SNP，跳过 indel、多等位位点和无 `GT` 位点；
- 不要求 VCF 已相位。

如果 header 没有 contig length，流程会使用该染色体最后一个保留 SNP 的位置作为扫描终点，因此推荐保留 contig length。

### 4.2 样本列表

`cultivar.list` 和 `wild.list` 均为每行一个 VCF 样本 ID：

```text
sample_001
sample_002
sample_003
```

空行和以 `#` 开头的行会被忽略。重复 ID 会去重。cultivar 与 wild 列表不能有交集，所有 ID 必须存在于 VCF header。

`--wild-list` 完全可选。不提供时自动进入 `CULTIVAR_ONLY`。

## 5. 最常用的运行命令

### 5.1 只有栽培种

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --outdir results_v17 \
  --prefix CA59 \
  2>&1 | tee CA59.v17.log
```

### 5.2 同时包含野生种和栽培种

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --wild-list wild.list \
  --outdir results_v17 \
  --prefix CA59 \
  2>&1 | tee CA59.v17.log
```

wild 不参与 block 边界或 cultivar 分群，只在最后计算到各 cultivar haplotype 的距离。

### 5.3 先测试一条染色体

```bash
python3 hapblock_v17.py \
  --vcf population.unfiltered.vcf.gz \
  --cultivar-list cultivar.list \
  --wild-list wild.list \
  --keep-chr Chr01 \
  --outdir test_Chr01 \
  --prefix CA59_Chr01
```

多条染色体可写为 `--keep-chr Chr01,Chr02,Chr03`。

## 6. 算法步骤

### 6.1 一次性基因型缓存

第一次运行时，V17 顺序读取一次 VCF，保存每条染色体的 SNP 坐标、`int8` 基因型、cultivar allele count 和 MAF。缓存默认位于：

```text
OUTDIR/01_genotype_cache/
```

同时生成 `OUTDIR/00_run_info/sample_qc.tsv`，记录每个样本的 call rate 和杂合率。缓存 manifest 保存 VCF 路径、文件大小、修改时间和样本列表；`--resume` 复用缓存前会检查这些信息。

### 6.2 Common-SNP backbone

默认使用 cultivar 中 `MAF >= 0.05` 和 `MAF >= 0.10` 的两个视图。高频变异更适合描述较古老、较稳定的群体结构。两个 MAF 视图可为同一边界提供独立支持。

### 6.3 多尺度窗口

默认 coarse 扫描：

| 窗口 | 步长 | 用途 |
|---:|---:|---|
| 1 Mb | 500 kb | 寻找大父 block |
| 500 kb | 250 kb | 支持父 block 和边界 |

默认 fine 扫描：

| 窗口 | 步长 | 用途 |
|---:|---:|---|
| 100 kb | 50 kb | 子 block |
| 50 kb | 25 kb | 子 block |
| 25 kb | 10 kb | 精细边界 |
| 10 kb | 10 kb | 最终边界分辨率 |

10/25 kb 窗口默认使用周围 100 kb 的上下文计算距离；有效 SNP 不足时扩大到 200 kb。这样可以报告 10 kb 级边界，而不会因为 10 kb 内 SNP 太少就产生大量假断点。

### 6.4 距离和分群

主距离为缺失校正后的 SNP/Mb：

```text
distance = mismatch_count / compared_site_count
           × eligible_site_count / block_length_Mb
```

默认 backbone 阈值为 `100 SNP/Mb`。完全连接聚类在该阈值切树，因此正式主要组内任意样本对的最大距离都不能超过阈值。聚类以后还会直接复核最大组内距离。

`--genotype-model inbred` 会把杂合型视为缺失，适合高度纯合材料；`--genotype-model dosage` 保留 0/1/2 状态，适合未相位二倍体。杂合较多的材料应使用 `dosage` 并仔细检查 PCA 和 callability。

### 6.5 父 block、短 gap 和子 block

相邻窗口只有在成员关系稳定时才合并，默认要求 ARI ≥ 0.85。低信息窗口可以桥接，但 coarse gap 同时受以下两个上限控制：

- `--coarse-max-gap-bp`，默认 1 Mb；
- `--max-gap-fraction`，默认当前候选 block 长度的 5%。

被桥接的区段写入 `*.internal_gaps.tsv`。持续的成员变化会切断 block。父 block 完成后，在父 block 内运行 fine 窗口并保留子 block；父 block 不会因存在子结构而被删除。

### 6.6 整个 block 的复核

窗口只负责提出候选边界。每个正式 block 都会在完整区间上重新抽取 SNP、计算距离和聚类。默认要求 coarse 候选至少由两个窗口尺度支持。这个步骤防止仅凭局部窗口产生不稳定 block。

### 6.7 RecentSubHap

在每个主要 backbone haplotype 内重新计算组内 MAF，并使用所有组内多态 SNP 检测近期结构。默认 recent 距离阈值为 `20 SNP/Mb`，比 backbone 更严格。

- 达到 `--min-recent-subhap-frequency` 和最小样本数的组标为 `SubHapXX`；
- 小组标为 `MinorSubHapXX` 并保留样本；
- SNP 信息不足时标为 `UnresolvedSubHap`，不纳入 recent-subhap 关联候选。

### 6.8 PCA

PCA 只用于 QC 和可视化。流程对 PCA 副本进行均值填补和位点标准化，并按最终 haplotype 着色。PCA 不覆盖完整连接聚类的正式样本归属。

## 7. 默认参数及其含义

| 参数 | 默认值 | 含义 |
|---|---:|---|
| `--backbone-mafs` | `0.05,0.10` | 定义 backbone 的 cultivar MAF 下限 |
| `--coarse-windows` | `1000000,500000` | 父 block 窗口 |
| `--fine-windows` | `100000,50000,25000,10000` | 子 block/边界窗口 |
| `--backbone-threshold-snp-per-mb` | `100` | backbone 最大组内距离 |
| `--recent-threshold-snp-per-mb` | `20` | recent subhap 最大组内距离 |
| `--membership-ari` | `0.85` | 相邻窗口成员一致性 |
| `--min-scale-support` | `2` | coarse block 最少支持尺度数 |
| `--min-structural-hap-frequency` | `0.05` | 正式结构 haplotype 最低频率 |
| `--min-structural-hap-samples` | `3` | 正式结构 haplotype 最少样本数 |
| `--min-association-hap-frequency` | `0.10` | 关联分析主要 haplotype 最低频率 |
| `--min-association-assigned-fraction` | `0.50` | 主要组合计最少覆盖样本比例 |
| `--min-association-haplotypes` | `2` | 关联 block 至少有两个主要组 |
| `--max-window-snps` | `2000` | 每个窗口距离计算的确定性 SNP 上限 |
| `--max-block-snps` | `10000` | 整个 block 复核的确定性 SNP 上限 |

SNP 超过上限时使用按坐标均匀、确定性的抽样，相同输入和参数会得到相同结果。

## 8. 针对不同目标的参数建议

### 8.1 几十 Mb 的着丝粒大 block

先使用默认 coarse 窗口。若大 block 被许多很短的低信息区打断，可逐步尝试：

```bash
--coarse-max-gap-bp 2000000 --max-gap-fraction 0.08
```

不要一次把 gap 放得很大。先查看 `*.internal_gaps.tsv`、`LOCAL_PASS_FRACTION`、PCA 和相邻 block 的样本成员是否确实一致。`max-gap-fraction` 仍能防止短 block 跨越过大的异常区。

### 8.2 更强调 10 kb 尺度

默认已经包含 10 kb 扫描。SNP 密度较高时可以缩小上下文：

```bash
--fine-context-bp 50000 --fine-max-context-bp 100000
```

若每 10 kb SNP 很少，应保留默认 100/200 kb 上下文，而不要降低 `--min-window-snps` 到极低值。

### 8.3 高纯合自交群体

使用默认：

```bash
--genotype-model inbred
```

先检查 `sample_qc.tsv`。若某些样本杂合率异常高，建议先核实样本、比对和 genotype calling，而不是直接放宽距离阈值。

### 8.4 杂合较多的二倍体群体

使用：

```bash
--genotype-model dosage
```

此模式比较 0/1/2 genotype 状态，但不会把 0|1 拆成两条单倍型。真正需要染色体级单倍型时，应先使用专门 phasing 工具并开发/使用 phased-haplotype 模式；V17.0.0 尚未实现该模式。

### 8.5 调整关联分析标准

用户提出的默认规则已经实现：每个主要 haplotype ≥10%，主要组合计覆盖 ≥50%，并且至少两个主要组。

例如要求更严格的 15%/70%：

```bash
--min-association-hap-frequency 0.15 \
--min-association-assigned-fraction 0.70
```

这些参数只决定哪些 block 进入关联矩阵，不会删除完整频率表中的小 haplotype。

## 9. 主要输出文件

### 9.1 Block 与层级

- `PREFIX.blocks.tsv`：所有 coarse 和 fine block 的主表，用 `BLOCK_LEVEL` 区分。
- `PREFIX.blocks.bed`：BED 形式的区间。
- `PREFIX.block_hierarchy.tsv`：fine block 到 coarse parent 的关系。
- `PREFIX.internal_gaps.tsv`：被桥接的内部低信息/异常区。
- `PREFIX.window_calls.tsv`：每个尺度、MAF 和窗口的详细状态。

`START0` 为 0-based，`END` 为右开区间终点，符合 BED。人类阅读区间为 `START0 + 1` 到 `END`。

`PREFIX.blocks.tsv` 中最重要的列：

- `N_MAJOR_HAPLOTYPES`：结构分析中达到 5%/最少样本条件的 haplotype 数；
- `ASSIGNED_FRACTION`：被分到这些结构主要组的 cultivar 比例；
- `MAX_WITHIN_DISTANCE_SNP_PER_MB`：主要组中最远样本对；
- `MIN_BETWEEN_DISTANCE_SNP_PER_MB`：不同主要组之间最近样本对；
- `BETWEEN_WITHIN_RATIO`：组间最小距离/组内最大距离；越大通常分离越清楚；
- `PAIR_CALLABLE_FRACTION`：样本对具有可比较位点的比例；
- `LOCAL_PASS_FRACTION`：构成该候选的局部窗口通过比例；
- `ASSOCIATION_READY`：是否达到关联分析 10%/50% 标准。

### 9.2 Haplotype 与样本

- `PREFIX.block_haplotype_frequency.tsv`：每个 block 的 haplotype 频率。
- `PREFIX.block_haplotype_samples.tsv`：长表，逐个样本列出 block 和 haplotype。
- `PREFIX.sample_by_block_haplotype.tsv`：宽表，每行一个样本、每列一个 block。

`HapXX` 是通过结构最低频率的组；`MinorXX` 是较小组。所有样本都保留，不会因是 minor 而从完整样本表删除。

### 9.3 Recent subhaplotype

- `PREFIX.backbone_to_recent_subhap.tsv`
- `PREFIX.recent_subhaplotype_frequency.tsv`
- `PREFIX.sample_by_recent_subhaplotype.tsv`
- `PREFIX.recent_subhap_association_status.tsv`
- `PREFIX.association_recent_subhap_onehot_matrix.tsv`

### 9.4 性状关联准备

- `PREFIX.association_ready_blocks.tsv`：通过默认或自定义规则的 block。
- `PREFIX.association_excluded_blocks.tsv`：未通过的 block及原因。
- `PREFIX.association_block_factor_matrix.tsv`：分类变量矩阵，适合把 block haplotype 当作多水平因子。
- `PREFIX.association_haplotype_onehot_matrix.tsv`：每个主要 haplotype 一列的 0/1 矩阵。

推荐优先使用 factor 模型对一个 block 做整体检验，再对显著 block 检查具体 haplotype。关联模型还应加入群体结构、亲缘关系和试验环境等协变量；V17 输出 haplotype 自变量，但不替代完整 GWAS/混合模型。

### 9.5 PCA、图形和报告

- `PREFIX.block_pca_scores.tsv`
- `PREFIX.block_pca_variance.tsv`
- `PREFIX.block_pca_top_loadings.tsv`
- `PREFIX.block_pca_plots.pdf`
- `PREFIX.block_landscape.pdf/png`
- `PREFIX.report.md`
- `PREFIX.output_index.tsv`

### 9.6 Wild 输出

仅在提供 `--wild-list` 时生成：

- `PREFIX.wild_to_cultivar_haplotype.tsv`：每个 wild 到各 cultivar haplotype 的距离及最近组。
- `PREFIX.wild_support_summary.tsv`：每个 cultivar haplotype 获得多少 wild 最近邻支持。

“wild 最近”只能解释为 ancestral-compatible 线索，不等于直接祖先证明。V17.0.0 不自动输出 ancestral/derived allele 极化结果。

## 10. 如何把结果与性状合并

表型表第一列应使用相同样本 ID。可将其与：

```text
PREFIX.association_block_factor_matrix.tsv
```

按 `SAMPLE` 合并。每个 block 是一个分类变量。若软件只接受数值自变量，则使用 one-hot 文件。分析前应确认：

1. 表型样本和 VCF 样本的交集；
2. 每个 haplotype 的实际样本数；
3. 群体结构是否与 haplotype 高度相关；
4. 是否需要只分析 coarse 或只分析 fine block，以控制多重检验；
5. parent 和 child 高度相关时不要把所有变量无条件同时放入同一普通线性模型。

## 11. 运行速度和资源

主要计算量约随样本数平方增长。缓存占用大致为：

```text
保留 SNP 数 × 所选样本数 × 1 byte
```

另加坐标和 MAF 数组。提高速度的顺序建议：

1. 先用 `--keep-chr Chr01` 检查参数；
2. 正式运行第一次完成缓存；同一 VCF 和列表再次运行时使用 `--resume`；
3. 只想检查 block 时加 `--skip-pca` 和/或 `--skip-recent-subhaps`；
4. 样本很多时先把 `--max-window-snps` 调为 1000、`--max-block-snps` 调为 5000 做参数探索；
5. 确认参数后用更高 SNP 上限复核最终结果；
6. 不要同时启动很多相同数据的大任务，以免磁盘缓存相互竞争。

复用缓存：

```bash
python3 hapblock_v17.py [原来的参数] --resume
```

如果 VCF 或样本列表变化，缓存校验会终止。使用新输出目录最安全；确实要覆盖时可显式加 `--force`。`--force` 会替换指定输出目录/缓存，请确认路径。

## 12. 如何判断结果是否可信

建议先检查：

1. `00_run_info/sample_qc.tsv`：低 call rate 或异常高杂合样本；
2. `blocks.tsv`：组内最大距离是否低于阈值、组间距离是否明显更大；
3. `window_calls.tsv`：边界两侧是否为持续变化，而不是单个低信息窗口；
4. `internal_gaps.tsv`：长 block 是否跨过过大的异常区；
5. PCA 图：最终 haplotype 是否具有可解释的分离，是否由单个离群样本驱动；
6. block 频率表：主要组样本数是否足够；
7. parent/child：显著 child 是否只是 parent 的重复编码。

第一次本地运行后，建议把以下文件发回用于一起判断准确性：

- `PREFIX.report.md`
- `PREFIX.blocks.tsv`
- `PREFIX.block_haplotype_frequency.tsv`
- `PREFIX.internal_gaps.tsv`
- `PREFIX.window_calls.tsv`
- `00_run_info/sample_qc.tsv`
- `PREFIX.block_landscape.png`
- 有 wild 时再加两个 wild 输出表。

## 13. 常见错误

### `Requested samples absent from VCF`

样本列表中的 ID 与 VCF header 不完全一致。检查大小写、不可见字符和重命名。

### `VCF is not chromosome-sorted`

同一染色体分成多个不连续段，或坐标倒序。先使用 bcftools 排序。

### block 很少

依次检查 common SNP 数、样本 call rate、`--min-window-snps`、距离阈值和 `--min-scale-support`。不要首先把阈值大幅放宽；先在一条染色体上查看 `window_calls.tsv` 的失败状态。

### block 太碎

检查是否由低 SNP 密度造成。可以适度增大 fine context、gap 上限或 ARI 容忍度。若成员确实持续改变，则碎分可能是真实子结构。

### block 过长

查看 `internal_gaps.tsv` 和 `LOCAL_PASS_FRACTION`。减小 `--coarse-max-gap-bp` 或 `--max-gap-fraction`，提高 `--membership-ari`，并确认两侧成员关系。

### 内存过高

优先减少同时分析的样本或按染色体运行；降低 `max-window-snps/max-block-snps` 主要减少计算时间和部分临时内存，但样本两两距离矩阵仍随样本数平方增长。

## 14. 软件测试与致谢

运行内置合成测试：

```bash
bash tests/run_tests.sh
```

测试覆盖完全连接避免链式合并、短 gap 桥接、持续边界切分、10%/50%关联规则、cultivar-only 模式和 cultivar+wild 模式。测试只使用合成数据，不会运行用户的大型 VCF。

HapBlock V17 的需求整理、流程重构、代码实现、测试和手册编写得到了 OpenAI Codex 的协助。最终生物学解释仍应由研究者结合物种遗传背景、样本设计和独立证据完成。
