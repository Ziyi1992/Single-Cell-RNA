# ============================================================
# LIVER SINGLE-CELL RNA-SEQ ANALYSIS WORKFLOW
# 肝脏单细胞RNA测序分析流程
#
# A profile-driven, modular workflow based on Seurat v5 and
# Bioconductor/OSCA principles
# 基于Seurat v5与Bioconductor/OSCA原则的配置驱动式模块化流程
# ============================================================
#
# PURPOSE / 目的
#
# This script provides a transparent and reproducible workflow for processing,
# quality assessment, cell-state discovery and sample-level statistical analysis
# of liver single-cell RNA-seq data. It is intended both for routine analysis and
# for methodological review by supervisors or collaborators.
#
# 本脚本用于对肝脏单细胞RNA测序数据进行透明、可重复的预处理、质量评估、
# 细胞状态识别及样本层面的统计分析。脚本既适用于日常数据分析，也便于导师或
# 合作者审阅每一步的分析依据、输入输出及质量控制标准。
#
# ANALYTICAL ROUTE / 分析路径
#
#   input → QC → doublet removal → normalization → HVG
#   → PCA/clustering/UMAP → cell annotation
#   → cell composition → pseudobulk differential expression
#   → functional interpretation/pathway activity → lineage subclustering/pseudotime
#   → exploratory cell-cell communication
#
#   数据输入 → 质量控制 → 双细胞识别与去除 → 标准化 → 高变基因筛选
#   → PCA/聚类/UMAP → 细胞类型注释
#   → 细胞组成分析 → pseudobulk差异表达分析
#   → 功能富集解释/通路活性分析 → 谱系内亚群/拟时分析
#   → 探索性细胞互作分析
#
# DIVISION OF METHODS / 方法分工
#
# Seurat v5 is used for object management, assay/layer handling, graph-based
# clustering, dimensionality reduction and visualization. Bioconductor/OSCA
# principles are used for sample-aware QC, scran normalization, variance-based
# feature selection and biological-replicate-level downstream inference.
#
# Seurat v5负责对象与assay/layer管理、图聚类、降维和可视化；
# Bioconductor/OSCA原则用于样本感知的质量控制、scran标准化、基于方差分解的
# 特征筛选，以及以生物学重复为单位的下游统计推断。
#
# STATISTICAL UNIT AND INFERENCE / 统计单位与推断原则
#
# Individual cells are observations for visualization and cell-state discovery,
# but they are not treated as independent biological replicates. Whenever group
# comparisons require formal inference, counts are aggregated by independent
# animal/patient and cell type before differential testing. This design avoids
# pseudoreplication caused by treating thousands of cells from the same specimen
# as thousands of independent experimental units.
#
# 单个细胞用于可视化和细胞状态识别，但不被视为彼此独立的生物学重复。
# 凡涉及组间正式统计推断，均先按独立动物/患者及细胞类型聚合计数，再进行
# 差异检验，从而避免将同一样本中的大量细胞误当作大量独立实验单位所造成的
# 伪重复（pseudoreplication）。
#
# INTERPRETATION BOUNDARY / 结果解读边界
#
# Clusters and UMAP coordinates are exploratory representations whose appearance
# depends on feature selection and analysis parameters. Cell-type labels must be
# supported by marker combinations, tissue context and biological knowledge;
# they should not be assigned from a single marker alone. Optional integration
# is used for sensitivity assessment and visualization, while raw RNA counts are
# retained for quantitative differential-expression analysis.
#
# 聚类和UMAP属于探索性表征，其结果会受到特征筛选和分析参数影响。细胞类型
# 注释必须综合marker组合、组织背景及生物学知识，不应依据单一marker下结论。
# 可选整合主要用于敏感性评估和可视化；定量差异表达分析始终保留并使用原始
# RNA计数，避免在整合表达矩阵上进行不恰当的统计推断。
#
# REPRODUCIBILITY AND CHECKPOINTS / 可重复性与检查点
#
# The workflow is divided into numbered modules. Each completed module writes
# tables, figures and logs to a dedicated directory and saves an RDS checkpoint
# when appropriate. A failed or revised module can therefore be rerun from the
# nearest checkpoint without repeating all upstream computation. The configured
# random seed, software log and input audit provide an explicit analysis record.
#
# 流程按编号模块组织。每个模块将表格、图片和日志写入独立目录，并在适当位置
# 保存RDS检查点。因此，某一步失败或需要调整时，可从最近的检查点继续运行，
# 无需重复全部上游计算。固定随机种子、软件环境日志和输入审计共同构成可追溯的
# 分析记录。
#
# ============================================================


# ============================================================
# REFERENCE FRAMEWORK / 参考框架
# ============================================================
#
# This template follows the analytical principles described in the
# Bioconductor Orchestrating Single-Cell Analysis (OSCA) workflow. The references
# below document the rationale for the core operations; the present script adapts
# those principles to a modular Seurat v5 workflow rather than reproducing one
# package-specific tutorial verbatim.
#
# 本模板遵循Bioconductor《Orchestrating Single-Cell Analysis》(OSCA)工作流的
# 分析原则。以下资料用于说明核心步骤的方法依据；本脚本将这些原则整合进模块化
# Seurat v5流程，并非逐行照搬某个软件包教程。
#
# Quality control / 质量控制
# https://bioconductor.org/books/release/OSCA.basic/quality-control.html
#
# Normalization / 标准化
# https://bioconductor.org/books/release/OSCA.basic/normalization.html
#
# Feature selection / 高变基因筛选
# https://bioconductor.org/books/release/OSCA.basic/feature-selection.html
#
# Dimensionality reduction / 降维
# https://bioconductor.org/books/release/OSCA.basic/dimensionality-reduction.html
#
# The script uses OSCA principles rather than copying one package workflow:
# 本脚本借鉴OSCA的分析原则，而不是机械照搬：
#
#   QC:
#     per-cell metrics + adaptive/sample-aware outlier detection
#     per-cell指标 + 自适应/sample-aware异常值判断
#
#   Normalization:
#     scran deconvolution size factors when OSCA_SCRAN is selected
#     选择OSCA_SCRAN时使用scran deconvolution normalization
#
#   Feature selection:
#     variance modelling and biological HVGs
#     variance modelling筛选生物学高变基因
#
#   Dimensionality reduction:
#     PCA on selected HVGs
#     在HVG上进行PCA
#
#   Seurat v5:
#     object management, graph clustering, UMAP and annotation
#     对象管理、图聚类、UMAP和人工注释
#
#   Statistical inference:
#     biological-replicate-level pseudobulk
#     以独立动物/患者为单位进行pseudobulk统计推断
#
# Three rules are applied throughout the workflow:
# 本流程始终遵循以下三项原则：
#
#   1. Preserve raw counts for count-based modelling.
#      保留原始计数，用于基于计数分布的统计模型。
#
#   2. Record QC decisions rather than silently discarding cells.
#      明确记录质量控制判定，不在无记录的情况下删除细胞。
#
#   3. Separate exploratory cell-level visualization from replicate-level
#      confirmatory inference.
#      区分细胞层面的探索性可视化与生物学重复层面的验证性统计推断。
#
# ============================================================


# ============================================================
# 0. USER SETTINGS / 用户设置
#
# This section is the user-editable analysis profile. It defines the dataset,
# experimental design, input route, module range and numerical parameters without
# altering the functions that implement the workflow.
#
# 本节是用户可编辑的分析配置区，用于定义数据集、实验设计、输入方式、运行模块
# 范围和数值参数，不需要修改后续实现分析流程的函数。
#
# Usually edit only this section and sample_metadata.csv.
# 通常更换数据集时只修改这一节和sample_metadata.csv。
#
# Before running a new dataset, verify every field against the experimental
# record. In particular, sample identity, biological replicate and condition must
# not be inferred from filenames unless that mapping has been independently
# confirmed.
#
# 更换数据集前，应依据实验记录逐项核对本节。尤其是样本身份、生物学重复和实验
# 分组，不应仅凭文件名推断，除非该对应关系已经得到独立确认。
# ============================================================

# Unique analysis identifier used in output folders, checkpoint names and logs.
# Use a stable identifier without path separators; changing it creates a separate
# result tree and does not rename existing outputs.
# 用于结果目录、检查点文件和日志的唯一分析标识。建议使用不含路径分隔符的稳定
# 名称；修改该值会建立新的结果目录，不会自动重命名已有结果。
DATASET_ID <- "GSE244475"

# ============================================================
# INPUT SETTINGS / 输入设置
# ============================================================
#
# Supported formats / 支持的输入格式：
#   "10X_MTX"    = matrix.mtx + features.tsv + barcodes.tsv
#   "10X_H5"     = 10x HDF5 file
#   "SEURAT_RDS" = prepared Seurat object
#   "SCE_RDS"    = SingleCellExperiment object
#
# Choose exactly one format. For RDS input, the object must retain a raw RNA
# count layer/matrix; normalized or integrated expression alone is insufficient
# for QC and pseudobulk differential analysis.
#
# 每次只选择一种格式。若读取RDS，对象必须保留RNA原始计数layer/矩阵；仅有
# 标准化或整合后的表达值不足以完成质量控制和pseudobulk差异分析。
# ============================================================
# 10X_MTX / 10X_H5 / SEURAT_RDS / SCE_RDS
INPUT_FORMAT <- "SEURAT_RDS"

# Directory containing the selected local input. Relative paths are resolved from
# the current R project directory; absolute paths are accepted for external data.
# 所选本地输入文件所在目录。相对路径以当前R项目目录为基准；外部数据也可使用
# 绝对路径。
INPUT_DIR <- "/Users/wangziyi/Desktop/R code/Proteomics&Transcriptomics/Single cell/GSE244475/GSE244475_standardized"

# TRUE reads files already present on disk. FALSE enables the GEO supplementary
# download route and interprets DATASET_ID as a GEO accession.
# TRUE表示读取磁盘上的本地文件；FALSE表示启用GEO补充文件下载流程，并将
# DATASET_ID解释为GEO登录号。
USE_LOCAL_FILES <- TRUE

# File used when INPUT_FORMAT = "SEURAT_RDS" or "SCE_RDS". The explicit
# file.path() construction avoids dependence on the current operating system's
# path separator.
# 当INPUT_FORMAT为"SEURAT_RDS"或"SCE_RDS"时使用的文件。采用file.path()
# 构建路径，可避免依赖特定操作系统的路径分隔符。
INPUT_RDS_FILE <- file.path(
  INPUT_DIR,
  "GSE244475_Liver_preQC.rds"
)

# File used only when INPUT_FORMAT = "10X_H5". Leave empty for other formats.
# 仅在INPUT_FORMAT = "10X_H5"时使用；其他输入格式可留空。
INPUT_H5_FILE <- ""

# ============================================================
# EXPERIMENTAL DESIGN / 实验设计
# ============================================================

# Exact condition labels expected in the metadata. CASE_GROUP is contrasted
# against CONTROL_GROUP in downstream models; therefore a positive log-fold
# change represents higher expression in CASE_GROUP unless a module explicitly
# states otherwise.
# metadata中预期出现的精确分组标签。下游模型以CASE_GROUP相对于CONTROL_GROUP
# 建立比较；除非具体模块另有说明，正向log-fold change表示CASE_GROUP中表达更高。
CONTROL_GROUP <- "Control"
CASE_GROUP <- "STZ"

# Species controls gene-symbol patterns and marker references. Mouse and human
# are supported through the complete workflow; rat is currently supported only
# through M10 because downstream annotation references require species-specific
# validation.
# 物种设置决定基因符号匹配规则和marker参考。mouse与human支持完整流程；rat目前
# 仅支持至M10，因为后续注释参考仍需进行物种特异性验证。
SPECIES <- "mouse"                 # mouse / human; rat is supported through M10 only
                                   # mouse / human；rat目前仅支持至M10

# Canonical tissue label used when optional tissue filtering is enabled. Its
# spelling must match the tissue value produced by the selected assignment mode.
# 启用可选组织过滤时使用的标准组织标签；其拼写必须与所选组织来源判定方式生成的
# tissue值完全一致。
TARGET_TISSUE <- "Liver"

# ============================================================
# MULTIPLEXING / 混样
# ============================================================

# Multiplexing strategy used to recover sample or tissue identity:
#   NONE                   = no demultiplexing is required;
#   HTO                    = classify cells from hashtag-oligo counts;
#   EXTERNAL_CELL_METADATA = import an externally validated cell-level mapping.
# 用于恢复样本或组织身份的混样策略：NONE表示无需拆分；HTO表示依据hashtag-oligo
# 计数分类；EXTERNAL_CELL_METADATA表示导入已经验证的细胞级对应表。
MULTIPLEXING_MODE <- "NONE"        # NONE / HTO / EXTERNAL_CELL_METADATA

# Biological meaning encoded by HTO labels. This determines which metadata field
# may be populated after demultiplexing; it does not infer label meaning from HTO
# names. Use SAMPLE or DONOR only when the experimental design supports that map.
# HTO标签所编码的生物学含义。该设置决定拆分后可填写哪个metadata字段，但不会
# 根据HTO名称自动猜测其含义。仅在实验设计明确支持时选择SAMPLE或DONOR。
HTO_LABEL_ROLE <- "TISSUE"          # TISSUE / SAMPLE / DONOR / OTHER

# ============================================================
# TISSUE ASSIGNMENT / 组织来源
# ============================================================

# Source used to assign tissue identity. SAMPLE_METADATA uses sample-level
# records; HTO uses demultiplexed labels; EXTERNAL_CELL_METADATA uses an imported
# cell-level table; NONE keeps existing tissue metadata without deriving it.
# 组织身份的判定来源。SAMPLE_METADATA使用样本级记录；HTO使用拆分标签；
# EXTERNAL_CELL_METADATA使用导入的细胞级表格；NONE不重新推导组织身份，仅保留
# 对象中已有的tissue信息。
TISSUE_ASSIGNMENT_MODE <- "NONE"
# SAMPLE_METADATA / HTO / EXTERNAL_CELL_METADATA / NONE

# If TRUE, retain only cells whose assigned tissue equals TARGET_TISSUE. Keep
# FALSE until tissue assignment has been audited, because an incorrect label map
# could otherwise remove valid cells irreversibly from the downstream object.
# 若为TRUE，仅保留组织标签等于TARGET_TISSUE的细胞。建议在组织来源映射完成审计
# 前保持FALSE，否则错误的标签对应关系可能使有效细胞从下游对象中被排除。
FILTER_TARGET_TISSUE <- FALSE

# ============================================================
# HTO SETTINGS / HTO设置
#
# These settings are retained for datasets that require HTO
# demultiplexing, but are not used for the current prepared RDS.
# 当前数据不会使用这些参数。
# ============================================================

# Accepted feature-type labels for candidate HTO assays in 10x data. Both names
# are retained because Cell Ranger exports may use either convention.
# 10x数据中可被识别为候选HTO assay的feature-type标签。保留两种名称，是因为
# 不同Cell Ranger输出可能采用不同命名方式。
HTO_FEATURE_TYPES <- c(
  "Antibody Capture",
  "Multiplexing Capture"
)

# Safety control: FALSE prevents an arbitrary non-gene-expression assay from
# being treated as HTO without explicit confirmation.
# 安全开关：FALSE可防止在未经明确确认的情况下，将任意非基因表达assay误识别为
# HTO数据。
ALLOW_NON_GEX_AS_HTO <- FALSE

# Explicit mapping from technical HTO names to biological labels. An empty vector
# means no remapping. Never complete this map by guessing from barcode order.
# HTO技术名称到生物学标签的明确对应表。空向量表示不进行重命名；不得依据barcode
# 顺序猜测并填写该对应关系。
# Example syntax only / 仅作语法示例：
# HTO_LABEL_MAP <- c("HTO1"="Liver", "HTO2"="Heart")
# 若HTO只是编码，必须依据真实metadata填写，代码不会猜。
HTO_LABEL_MAP <- character(0)

# Positive quantile used by HTO demultiplexing. Higher values are more stringent
# and may increase the number of negative/unclassified cells; its suitability
# must be checked from HTO count distributions and classification summaries.
# HTO拆分使用的阳性分位数阈值。数值越高通常越严格，也可能增加negative或未分类
# 细胞；应结合HTO计数分布和分类汇总判断是否合适。
HTO_POSITIVE_QUANTILE <- 0.99

# TRUE retains singlets for downstream analysis after HTO classification and
# excludes HTO doublets/negatives. Classification totals are saved before removal.
# TRUE表示HTO分类后仅保留singlet进入下游，并排除HTO doublet/negative；删除前会
# 保存各类别数量，便于审计。
KEEP_HTO_SINGLETS_ONLY <- TRUE

# Path to a validated cell-level annotation table used only for
# EXTERNAL_CELL_METADATA mode. Cell identifiers must match object barcodes.
# 仅在EXTERNAL_CELL_METADATA模式下使用的、已经验证的细胞级注释表路径；其中的
# 细胞标识必须与对象barcode一致。
EXTERNAL_CELL_METADATA_FILE <- ""

# ============================================================
# METADATA COLUMNS / metadata列名
# ============================================================

# Column containing the sequencing/library sample identifier. This may describe
# a technical library and is not automatically equivalent to an independent
# biological replicate.
# 测序样本或文库标识所在列。该标识可能仅代表技术文库，不一定等同于独立的
# 生物学重复。
SAMPLE_ID_COL <- "sample_id"

# Column identifying the independent animal, patient or experimental specimen.
# This is the default unit for pseudobulk aggregation and is therefore essential
# for valid group-level inference.
# 独立动物、患者或实验标本的标识列。该列默认作为pseudobulk聚合单位，是保证
# 组间统计推断有效性的关键字段。
BIOLOGICAL_REPLICATE_COL <- "biological_replicate"

# Column containing experimental group labels. Values must include the exact
# CONTROL_GROUP and CASE_GROUP strings defined above.
# 实验分组标签所在列；其取值必须包含上方定义的CONTROL_GROUP和CASE_GROUP精确
# 字符串。
CONDITION_COL <- "condition"

# Column containing tissue identity. Required only when tissue assignment or
# target-tissue filtering is requested, but retained in all outputs when present.
# 组织来源所在列。仅在进行组织判定或目标组织过滤时强制需要；若原数据已包含，
# 则会保留在所有下游输出中。
TISSUE_COL <- "tissue"

# Optional sample-level metadata table. Leave empty to use the default file under
# config/. One row should represent one sample-level record, with unique keys and
# experimentally verified mappings.
# 可选的样本级metadata表路径。留空时使用config/目录下的默认文件。每行应代表一条
# 样本记录，键值必须唯一，且对应关系应由实验记录确认。
SAMPLE_METADATA_FILE <- ""

# ============================================================
# SOFTWARE / 软件
# ============================================================

# TRUE allows the environment module to install missing packages. For a managed
# server or a frozen reproducible environment, set FALSE and install dependencies
# separately so that package versions remain under external control.
# TRUE允许环境检查模块安装缺失的软件包。在受管理的服务器或固定版本环境中，建议
# 设为FALSE并单独安装依赖，以便由外部环境统一控制软件版本。
INSTALL_PACKAGES <- TRUE

# Parallel execution is optional and affects runtime rather than scientific logic.
# N_WORKERS should not exceed the compute resources allocated to the R session.
# 并行计算仅影响运行时间，不改变科学分析逻辑。N_WORKERS不应超过当前R会话实际
# 分配到的计算资源。
USE_PARALLEL <- FALSE
N_WORKERS <- 4L

# ============================================================
# MODULES / 模块
# ============================================================

# Inclusive module range for the current run. On a new dataset, first run M00-M03
# and review environment, metadata and input-audit reports before reading and
# filtering the full object. To resume, RUN_FROM should point to the next module
# and the required upstream checkpoint must already exist.
# 本次运行模块的闭区间。新数据集建议先运行M00-M03，并检查环境、metadata和输入
# 审计报告，再读取和过滤完整对象。续跑时，RUN_FROM应指向下一模块，且所需上游
# 检查点必须已经存在。
RUN_FROM <- 15L
RUN_TO <- 15L

# Optional integration is a sensitivity branch attached to M10. It does not
# replace the primary non-integrated object and is not used as the count source
# for pseudobulk differential-expression testing.
# 可选整合是连接在M10上的敏感性分析分支，不替代主分析的未整合对象，也不作为
# pseudobulk差异表达检验的计数来源。
RUN_OPTIONAL_INTEGRATION <- FALSE

# QC / 质量控制
# OSCA applies adaptive, sample-aware outlier detection; FIXED applies only the
# explicit numerical thresholds below; HYBRID combines both and is therefore more
# stringent. Whichever mode is selected, flags and summaries should be reviewed
# before interpreting the retained population.
# OSCA使用自适应、样本感知的异常值判定；FIXED仅使用下方固定阈值；HYBRID同时
# 使用两者，因此通常更严格。无论选择哪种模式，都应先审阅标记和汇总，再解释
# 最终保留的细胞群体。
QC_MODE <- "OSCA"                  # OSCA / FIXED / HYBRID

# TRUE removes cells that fail the selected QC rule. FALSE calculates and saves
# all QC flags without filtering, which is useful for threshold review.
# TRUE表示删除未通过所选QC规则的细胞；FALSE仅计算并保存所有QC标记而不实施过滤，
# 适合用于阈值审查。
APPLY_QC_FILTER <- TRUE

# Fixed thresholds used in FIXED and HYBRID modes. Inf disables the corresponding
# upper bound. Mitochondrial percentage is expressed on a 0-100 scale.
# FIXED和HYBRID模式使用的固定阈值。Inf表示不启用对应上限；线粒体比例采用0-100
# 百分数尺度。
QC_MIN_FEATURES <- 200
QC_MAX_FEATURES <- Inf
QC_MIN_COUNTS <- 500
QC_MAX_COUNTS <- Inf
QC_MAX_MT <- 30

# Metadata column defining QC batches for adaptive outlier detection. A sample
# identifier is usually appropriate because library quality distributions may
# differ across preparations. Each batch must contain enough cells for stable
# estimates.
# 自适应异常值判定所使用的QC批次列。通常选择样本标识，因为不同文库制备的质量
# 分布可能不同；每个批次需有足够细胞才能稳定估计。
QC_BATCH_COL <- SAMPLE_ID_COL

# TRUE preserves per-cell decision flags and reason summaries, allowing every
# exclusion to be traced after filtering.
# TRUE表示保留逐细胞QC判定及原因汇总，使过滤后仍可追溯每个排除决定。
SAVE_QC_FLAGS <- TRUE

# Doublets / 双细胞
# The first switch runs scDblFinder; the second controls whether predicted
# doublets are removed. Keeping detection and removal separate permits diagnostic
# review before exclusion.
# 第一个开关控制是否运行scDblFinder；第二个开关控制是否删除预测双细胞。将检测与
# 删除分开，便于在排除前先审阅诊断结果。
RUN_SCDOUBLETFINDER <- TRUE
REMOVE_DOUBLETS <- TRUE

# Maximum number of cells sampled from each sample-by-class combination for the
# M08 complexity diagnostic. The complete classification is retained in the
# object and summary table; subsampling affects visualization only and prevents
# large datasets from producing an unreadable overplotted panel.
# M08文库复杂度诊断图中，每个sample×classification组合最多抽取的细胞数。
# 完整分类仍保留在对象和汇总表中；抽样仅影响绘图，用于避免大数据集严重遮盖。
DOUBLET_DIAGNOSTIC_MAX_CELLS_PER_SAMPLE_CLASS <- 1500L

# Normalization / HVG
# OSCA_SCRAN uses deconvolution size factors and variance modelling; LOGNORMALIZE
# uses Seurat log-normalization; SCT uses SCTransform. The selected branch changes
# preprocessing but must preserve raw RNA counts for later pseudobulk analysis.
# OSCA_SCRAN使用deconvolution size factor和方差建模；LOGNORMALIZE使用Seurat
# 对数标准化；SCT使用SCTransform。所选分支会改变预处理方式，但必须保留RNA原始
# 计数以供后续pseudobulk分析。
NORMALIZATION_METHOD <- "OSCA_SCRAN"  # OSCA_SCRAN / LOGNORMALIZE / SCT

# Target number of highly variable genes used for PCA. This controls the feature
# space for exploratory structure discovery, not the genes available for DE.
# 用于PCA的目标高变基因数量。该参数控制探索性结构识别的特征空间，不限制差异
# 表达分析可使用的基因范围。
N_HVG <- 3000L

# Number of leading HVGs labelled in the M09 diagnostic figures. Labels are
# selected from the fitted biological-variance component and are for display
# only; changing this value does not alter the HVGs used for PCA.
# M09诊断图中标注的前列HVG数量。标签依据拟合得到的biological variance component
# 选择，仅影响图片展示，不改变实际进入PCA的HVG集合。
HVG_DIAGNOSTIC_LABEL_TOP <- 12L

# If TRUE, mitochondrial percentage is regressed during scaling. Use only with a
# clear technical rationale, because mitochondrial signal may also reflect real
# biological state.
# 若为TRUE，在scale阶段回归线粒体比例。仅在具有明确技术依据时启用，因为线粒体
# 信号也可能反映真实的生物学状态。
REGRESS_PERCENT_MT <- FALSE

# ============================================================
# PCA / CLUSTERING / UMAP
# PCA / 聚类 / UMAP
# ============================================================

# Maximum number of principal components calculated. Computing more PCs than are
# used allows the elbow and variance-explained diagnostics to be inspected.
# 实际计算的主成分上限。计算数量高于正式使用数量，可用于检查elbow plot和方差
# 解释度诊断。
N_PCS_COMPUTE <- 50L

# Number of PCs used for neighbor graph and UMAP.
# This should be checked against the PCA elbow plot.
# 用于邻居图、聚类和UMAP的PC数量。
# 需要结合PCA elbow plot检查，不应机械套用。
N_PCS_USE <- 30L

# Number of nearest neighbors for graph construction.
# 构建邻居图时使用的近邻数量。
K_NEIGHBORS <- 20L

# Final clustering resolution used for downstream analysis.
# 下游正式分析采用的聚类resolution。
CLUSTER_RESOLUTION <- 0.5

# Explicit graph names avoid hard-coding assay-specific names
# such as "RNA_snn".
#
# 显式定义graph名称，避免代码依赖RNA_snn等
# Seurat自动生成的assay-specific名称。
CLUSTER_NN_GRAPH <- "analysis_nn"
CLUSTER_SNN_GRAPH <- "analysis_snn"

# Resolution grid is used ONLY for sensitivity diagnostics.
# 最终正式结果仍由CLUSTER_RESOLUTION决定。
CLUSTER_RESOLUTION_GRID <- c(
  0.2,
  0.4,
  0.6,
  0.8,
  1.0
)

UMAP_N_NEIGHBORS <- 30L
UMAP_MIN_DIST <- 0.3

# Optional integration / 可选整合
# Currently supported strategy for the optional Seurat integration branch. Its
# assumptions and output must be evaluated independently of the primary analysis.
# 可选Seurat整合分支所使用的策略；其假设和结果必须独立于主分析进行评估。
INTEGRATION_METHOD <- "CCA"

# Annotation / 注释
# TRUE creates a cluster-level manual-annotation template. The generated file is
# intentionally incomplete until cell_type is filled after reviewing markers,
# expression plots and biological context, then M11 is rerun to write labels back.
# TRUE表示生成cluster层面的人工注释模板。该文件会有意保持未完成状态，需在审阅
# marker、表达图和生物学背景后填写cell_type，再重新运行M11写回标签。
CREATE_MANUAL_ANNOTATION_TEMPLATE <- TRUE

# Maximum number of cells sampled from each displayed category for the
# cell-level ComplexHeatmap. Downsampling prevents abundant populations from
# dominating image width while retaining within-category heterogeneity.
# ComplexHeatmap中每个展示分类最多抽取的cell数量。下采样可以避免高丰度群体占据
# 绝大部分图片宽度，同时保留每个分类内部的表达异质性。
COMPLEX_HEATMAP_MAX_CELLS_PER_CATEGORY <- 150L

# Symmetric clipping limit for per-gene Z-scores in ComplexHeatmap. Values
# outside this range are clipped for visualization only.
# ComplexHeatmap中逐基因Z-score的对称截断范围；超出范围的值仅在绘图时截断。
COMPLEX_HEATMAP_Z_LIMIT <- 2

# Composition / 细胞组成
# Minimum cell count required for a sample-by-cell-type entry to appear in the
# composition summary. This is a reporting threshold, not a substitute for a
# replicate-level compositional model.
# 样本×细胞类型组合进入组成汇总所需的最少细胞数。该值是报告阈值，不能替代以
# 生物学重复为单位的组成统计模型。
MIN_CELLS_FOR_COMPOSITION <- 1L

# TRUE performs an exploratory replicate-level comparison of cell-type
# proportions between CONTROL_GROUP and CASE_GROUP. This screening summary is
# useful for visualization but is not a replacement for a dedicated
# compositional differential-abundance model.
# TRUE表示对CONTROL_GROUP和CASE_GROUP之间的细胞类型比例进行探索性的重复层面
# 比较。该结果适合筛查和可视化，但不能替代专门的组成型差异丰度模型。
RUN_EXPLORATORY_COMPOSITION_TEST <- TRUE

# Minimum number of biological replicates required in each condition before an
# exploratory proportion test is attempted for a cell type.
# 每个condition至少需要达到该数量的biological replicate，才对相应细胞类型执行
# 探索性比例检验。
COMPOSITION_TEST_MIN_REPLICATES_PER_GROUP <- 2L

# Multiple-testing correction applied across cell types in the exploratory
# composition screen. "BH" controls the false-discovery rate.
# 探索性组成筛查中跨细胞类型进行多重检验校正的方法；"BH"用于控制错误发现率。
COMPOSITION_P_ADJUST_METHOD <- "BH"

# Pseudobulk DE / Pseudobulk差异表达
# Pseudobulk aggregates raw counts within each biological replicate and cell type,
# then fits count-based group comparisons. The minimum-cell threshold protects
# against extremely sparse aggregates; the replicate threshold is checked per
# group and cell type before a model is fitted.
# Pseudobulk先在每个生物学重复和细胞类型内聚合原始计数，再建立基于计数的组间
# 比较。最少细胞阈值用于避免极度稀疏的聚合样本；建模前还会按组别和细胞类型
# 检查生物学重复数量。
RUN_PSEUDOBULK <- TRUE
MIN_CELLS_PER_PSEUDOBULK <- 20L
MIN_REPLICATES_PER_GROUP <- 2L
PSEUDOBULK_REPLICATE_COL <- BIOLOGICAL_REPLICATE_COL

# False-discovery-rate threshold used to classify and report significant genes
# after edgeR testing. This threshold affects result summarization and plotting;
# it does not alter model fitting or the complete edgeR result tables.
# edgeR检验完成后用于判定和汇报显著基因的FDR阈值。该阈值只影响结果汇总和绘图，
# 不改变模型拟合过程，也不会删除完整edgeR结果表中的其他基因。
PSEUDOBULK_DE_FDR_THRESHOLD <- 0.05

# Minimum absolute log2 fold change required for the significant-DEG summary.
# Positive values indicate higher expression in CASE_GROUP; negative values
# indicate higher expression in CONTROL_GROUP.
# 显著DEG汇总所要求的最小绝对log2 fold change。正值表示CASE_GROUP表达更高，
# 负值表示CONTROL_GROUP表达更高。
PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD <- 0.5

# Maximum number of significant genes labelled in each direction for every cell
# type. Labels are selected by FDR first and effect size second to keep combined
# figures readable without changing which genes are classified as significant.
# 每种细胞类型、每个变化方向最多标注的显著基因数量。标签优先依据FDR、其次依据
# 效应量选择，以控制合并图片的拥挤程度，但不会改变显著基因的判定结果。
PSEUDOBULK_DE_LABELS_PER_DIRECTION <- 3L

# Optional sample-level covariates included in the DE design, for example sex or
# batch. Every covariate must vary independently enough to be estimable and must
# not be perfectly confounded with condition.
# 可选的样本级差异模型协变量，例如sex或batch。每个协变量必须具有足够独立变异，
# 且不能与condition完全混杂，否则模型无法可靠估计。
DE_COVARIATES <- character(0)

# Functional interpretation of pseudobulk DE / Pseudobulk差异基因功能解释
#
# This branch reads the complete M13 edgeR tables; it never recomputes cell-level
# differential expression. Significant genes are classified by annotated cell
# type and direction, displayed in a cross-cell-type log2FC heatmap, and tested
# for GO Biological Process and KEGG pathway over-representation using the genes
# that passed edgeR filtering in the same cell type as the statistical universe.
# 本分支直接读取M13的完整edgeR结果，不会重新进行cell层面差异检验。
# 显著基因按已注釄cell type及变化方向分类，并通过跨cell type的log2FC
# 热图展示；GO Biological Process与KEGG过度富集使用同一cell type内通过
# edgeR表达过滤的全部基因作为检验背景，而不是使用整个注释数据库。
RUN_FUNCTIONAL_INTERPRETATION <- TRUE

# Number of source DEGs retained in each cell-type-by-direction heatmap block.
# Selection is based on FDR first and absolute log2FC second; it affects display
# only and never changes the complete DEG tables or enrichment input.
# 差异基因分类热图中每个cell type×变化方向最多展示的来源DEG数。
# 先按FDR、再按绝对log2FC选择；该参数只影响图片，不改变完整结果表及富集输入。
DEG_CLASS_HEATMAP_TOP_PER_CLASS <- 10L

# A cell-type/direction class needs at least this many successfully mapped
# significant genes before enrichment is attempted. Classes below the threshold
# are retained in an explicit audit table instead of failing the module.
# 每个cell type/方向类别至少需要达到该数量的成功映射显著基因，才执行富集分析。
# 低于门槛的类别会记录在审计表中，不会导致整个模块失败。
ENRICHMENT_MIN_MAPPED_DEGS <- 5L
ENRICHMENT_FDR_THRESHOLD <- 0.05
ENRICHMENT_TOP_TERMS_PER_CLASS <- 8L

# ggkegg is generated only for FDR-significant KEGG pathways. Selection is first
# performed within every cell-type-by-direction class so one highly significant
# class cannot monopolize all pathway panels. GGKEGG_MAX_PATHWAYS then provides a
# global safety cap for runtime and output volume. No pathway name or ID is fixed.
# ggkegg只对FDR显著的KEGG通路作图。通路首先在每个cell type×方向
# 类别内选取，避免某一高显著类别占用全部图片名额；随后再使用
# GGKEGG_MAX_PATHWAYS限制总运行时间及输出数量。不写死通路名称或ID。
RUN_GGKEGG <- TRUE
GGKEGG_PATHWAYS_PER_CLASS <- 2L
GGKEGG_MAX_PATHWAYS <- 12L

# Maximum mapped DEGs shown in the companion effect-size panel beside each KEGG
# map. All mapped genes remain in the output table; this limit affects labels and
# plotting only. Pathway maps preserve the original KEGG labels and connections.
# 每张KEGG通路图旁的效应量panel最多展示的映射DEG数。所有映射基因
# 仍完整写入表格；该参数只控制标签和图片密度。通路图保留KEGG原始
# 标签及连接结构。
GGKEGG_MAX_GENES_IN_COMPANION <- 15L

# Pathway activity / 通路活性
#
# GSVA is performed on TMM-normalized pseudobulk logCPM matrices retained from
# M13. Each biological replicate remains one column and each annotated cell type
# is analysed independently. The resulting pathway scores are compared with
# limma using the same condition and optional covariates as pseudobulk edgeR.
# GSVA使用M13保留的TMM标准化pseudobulk logCPM矩阵。每个生物学重复
# 始终是一列，每种已注释cell type分别分析。通路score使用limma按与
# pseudobulk edgeR一致的condition及可选协变量设计进行组间比较。
RUN_GSVA <- TRUE

# KEGG pathways are retained only when the number of expressed member genes in a
# cell type lies within this interval. Very small sets are unstable, whereas very
# broad sets are difficult to interpret and may dominate rank-based scores.
# 仅当某条KEGG通路在当前cell type中可检测成员基因数位于该区间时才进入
# GSVA。过小基因集的score不稳定，过宽基因集则较难进行具体生物学解释。
GSVA_MIN_GENE_SET_SIZE <- 10L
GSVA_MAX_GENE_SET_SIZE <- 500L

# Prior count used only when converting the filtered TMM-normalized pseudobulk
# counts to continuous logCPM values for GSVA. It does not alter the edgeR model.
# 将经过滤并完成TMM标准化的pseudobulk counts转换为GSVA所需连续
# logCPM时使用的prior count。该参数不会改变edgeR统计模型。
GSVA_LOGCPM_PRIOR_COUNT <- 2

# FDR and minimum absolute pathway-score difference used for the concise
# significant-pathway table. Complete limma results are always retained.
# 简化显著通路表所使用的FDR及最小绝对通路score差异。完整limma结果
# 始终保留，不受该汇报阈值删减。
GSVA_FDR_THRESHOLD <- 0.05
GSVA_MIN_ABS_SCORE_DIFFERENCE <- 0.10

# Number of pathways displayed per cell type in sample-score panels and heatmaps.
# Selection uses FDR first and absolute score difference second; it changes
# presentation only and does not affect pathway testing.
# 每个cell type的sample-level score图及热图最多展示的通路数。先按FDR、
# 再按绝对score差异选取；该设置仅影响图片，不影响通路检验。
GSVA_TOP_PATHWAYS_PER_CELL_TYPE <- 8L

# Maximum number of pathways in the compact cross-cell-type overview. The
# pathway-by-cell-type bubble matrix and effect heatmap use the same ranked set,
# preventing an excessively tall figure when many cell types are present.
# 跨cell type紧凑总览图最多展示的通路数。气泡矩阵与效应热图使用同一组
# 排序后的通路，避免在cell type较多时生成不适合论文正文的超长图片。
GSVA_OVERVIEW_MAX_PATHWAYS <- 30L

# Final evidence synthesis / 最终证据整合
#
# TRUE asks M15 to consolidate completed module outputs into an auditable final
# evidence matrix, statistically selected key-gene tables, targeted pathway-axis
# summaries and a deterministic bilingual report. This step reads existing result
# tables only: it does not refit models, change significance calls or use an AI
# model to generate biological conclusions.
# TRUE表示由M15将已完成模块的结果整合为可审计的最终证据矩阵、统计学筛选的关键
# 基因表、定向通路轴汇总及规则驱动的双语报告。本步骤只读取现有结果表，不会重新
# 拟合模型、改变显著性判定，也不会调用AI模型自动生成生物学结论。
RUN_FINAL_EVIDENCE_SYNTHESIS <- TRUE

# Maximum number of significant genes displayed for each cell type in the compact
# final figure. Complete statistically significant genes remain in the exported
# table; selection uses FDR first and absolute log2 fold change second.
# 最终紧凑图片中每种cell type最多展示的显著基因数。完整显著基因始终保存在输出
# 表格中；图片依次按FDR和绝对log2 fold change选择，避免人工挑选基因。
FINAL_EVIDENCE_TOP_GENES_PER_CELL_TYPE <- 8L

# Common adjusted-P-value threshold used only when M15 classifies evidence from
# modules that do not already provide a logical significance column. Module-
# specific thresholds remain authoritative for DEG, enrichment and GSVA tables.
# 仅当上游结果表未提供逻辑型显著性列时，M15使用该统一校正P值阈值进行证据分类；
# DEG、富集和GSVA仍以各自模块的专用阈值为准。
FINAL_EVIDENCE_FDR_THRESHOLD <- 0.05

# A nominal composition P value below this threshold, when the adjusted P value
# is not significant, is labelled only as an exploratory trend. It never counts
# as a statistically significant evidence family in the final synthesis.
# 当细胞比例校正后P值不显著、但原始P值低于该阈值时，仅标记为探索性趋势；该趋势
# 不会在最终整合中计作显著证据。
FINAL_EVIDENCE_COMPOSITION_TREND_P <- 0.05

# Optional hypothesis-led KEGG axes for the secondary focused summary. Pathway
# identifiers are database-defined KEGG core IDs and therefore remain compatible
# across supported species prefixes (for example mmu, hsa or rno). These axes do
# not influence primary DEG/pathway discovery or the ranking of key cell types.
# 用于二级定向汇总的可选KEGG研究轴。此处使用数据库定义的KEGG核心ID，因此可兼容
# 不同物种前缀（如mmu、hsa或rno）。这些研究轴不会影响全局DEG/通路发现，也不会
# 影响关键cell type的排序。
FINAL_EVIDENCE_TARGETED_KEGG_AXES <- list(
  glycogen_glucose = c(
    "00500", # Starch and sucrose metabolism / 淀粉和蔗糖代谢
    "00010", # Glycolysis / Gluconeogenesis / 糖酵解与糖异生
    "04910", # Insulin signaling / 胰岛素信号
    "04911", # Insulin secretion / 胰岛素分泌
    "04922", # Glucagon signaling / 胰高血糖素信号
    "04068"  # FoxO signaling / FoxO信号
  ),
  energy_sensing = c(
    "04931", # Insulin resistance / 胰岛素抵抗
    "04152", # AMPK signaling / AMPK信号
    "04150", # mTOR signaling / mTOR信号
    "00190"  # Oxidative phosphorylation / 氧化磷酸化
  ),
  lipid_MASLD = c(
    "01212", # Fatty acid metabolism / 脂肪酸代谢
    "00071", # Fatty acid degradation / 脂肪酸降解
    "00062", # Fatty acid elongation / 脂肪酸延长
    "01040", # Unsaturated fatty-acid biosynthesis / 不饱和脂肪酸合成
    "04932"  # Non-alcoholic fatty liver disease / 脂肪性肝病
  ),
  fibrosis_ECM = c(
    "04512", # ECM-receptor interaction / ECM-受体互作
    "04518", # Integrin signaling / Integrin信号
    "04610"  # Complement and coagulation cascades / 补体与凝血
  ),
  inflammation_checkpoint = c(
    "04657", # IL-17 signaling / IL-17信号
    "05235"  # PD-L1 and PD-1 checkpoint / PD-L1与PD-1检查点
  ),
  HCC_related = c(
    "05225", # Hepatocellular carcinoma / 肝细胞癌
    "05200", # Pathways in cancer / 癌症通路
    "05205"  # Proteoglycans in cancer / 癌症中的蛋白聚糖
  )
)

FINAL_EVIDENCE_TARGETED_AXIS_LABELS <- c(
  glycogen_glucose =
    "Glycogen and glucose metabolism",
  energy_sensing =
    "Energy sensing",
  lipid_MASLD =
    "Lipid metabolism and MASLD-related pathways",
  fibrosis_ECM =
    "Fibrosis and extracellular matrix",
  inflammation_checkpoint =
    "Inflammation and immune checkpoint",
  HCC_related =
    "HCC-related pathways"
)

# Within-lineage subcluster reporting / 谱系内亚群报告
# M14B already recomputes a graph within every configured lineage. These settings
# add an exploratory marker table, dot plot and heatmap for those expression-
# derived subclusters. Numeric subclusters are deliberately not auto-renamed as
# biological states; the marker audit supports subsequent manual interpretation.
# M14B已会在每条配置谱系内重新构建图并聚类。以下参数用于新增探索性
# 亚群marker表、气泡图和热图。数字亚群不会被自动命名为生物学状态；
# marker审计结果用于后续人工判读。
LINEAGE_SUBCLUSTER_MARKER_MIN_PCT <- 0.10
LINEAGE_SUBCLUSTER_MARKER_LOGFC <- 0.25
LINEAGE_SUBCLUSTER_TOP_MARKERS <- 5L

# Exploratory cell-cell communication / 探索性细胞互作
#
# CellChat infers candidate communication from normalized ligand/receptor
# co-expression and a curated prior database. It does not directly demonstrate
# physical contact, ligand secretion, receptor activation or causal signalling.
# Separate models are fitted for CONTROL_GROUP and CASE_GROUP; labels are read
# from the analysis profile and are never fixed to a particular experiment.
# CellChat基于标准化后的配体/受体共表达及经整理的先验数据库推断候选互作。
# 它不能直接证明物理接触、配体分泌、受体激活或因果信号。流程分别对
# CONTROL_GROUP与CASE_GROUP建模，实际组名从用户配置读取，不固定于某个实验。
RUN_CELLCHAT <- TRUE

# A cell type is compared only when it reaches both the cell-count and independent-
# replicate thresholds in BOTH conditions. Requiring a common cell-type set makes
# network differences interpretable and avoids a one-sample population defining
# an apparent condition-specific interaction.
# 只有在两个condition中都同时达到细胞数和独立生物学重复门槛的cell type，
# 才进入组间互作比较。两组使用共同cell type集，避免单一样本中的稀有群体
# 产生表面上的condition特异互作。
CELLCHAT_MIN_CELLS_PER_TYPE_CONDITION <- 30L
CELLCHAT_MIN_REPLICATES_PER_TYPE_CONDITION <- 2L

# CellChat's default triMean estimator is retained. min.cells is an additional
# within-model group-size filter; top-pair limits affect visualization only.
# 保留CellChat默认triMean表达估计。min.cells是模型内额外的群体大小过滤；
# top pair数量只影响图片展示，不删除完整互作结果表。
CELLCHAT_MIN_CELLS_FOR_INTERACTION <- 10L
CELLCHAT_TOP_LR_PAIRS <- 40L
CELLCHAT_NETWORK_TOP_EDGES <- 40L
CELLCHAT_DATABASE_CATEGORIES <- c(
  "Secreted Signaling",
  "ECM-Receptor",
  "Cell-Cell Contact"
)

# Optional focused analysis / 可选重点分析
# Optional subsets for targeted reporting after the global analysis. Empty vectors
# mean no predefined focus; they do not remove cells or genes from upstream steps.
# 全局分析完成后用于重点报告的可选细胞类型和基因集合。空向量表示不预设重点，
# 不会从上游步骤中删除任何细胞或基因。
FOCUS_CELL_TYPES <- character(0)
FOCUS_GENES <- character(0)

# Optional lineage-restricted pseudotime analysis / 可选的谱系限定拟时分析
#
# Pseudotime is suitable only for cells that may belong to one biologically
# plausible lineage or state-transition system. It must not be applied to all
# annotated cell types simply because they are present in the same Seurat object.
# The predefined profiles below are evaluated independently. A profile is skipped
# when none of its requested cell types are present, so enabling this module does
# not force unrelated cell types into one trajectory.
#
# 拟时分析只适用于可能属于同一生物学谱系或状态转换体系的细胞。不能仅因为不同
# cell type同时存在于一个Seurat对象中，就把全部细胞放入同一条轨迹。通用配置中
# 下列预设谱系彼此独立评估。当某个数据集不存在目标cell type时，该谱系会被跳过，
# 因此启用本模块不会把不相关的细胞强行连接到同一条轨迹。
RUN_PSEUDOTIME <- TRUE

# Five priority profiles relevant to hepatic glycogen biology are defined here.
# profile_id is used only for stable output-directory names. target_cell_types
# defines the cells entering one independent analysis; root/endpoint identities
# constrain orientation but do not prove differentiation. Feature genes are
# visualized after fitting and never determine the trajectory itself.
#
# 这里预设五条与肝糖原研究相关、彼此独立的优先谱系。profile_id仅用于生成稳定的
# 输出目录名；target_cell_types定义每次单独分析所纳入的细胞；root/endpoint身份
# 只约束方向，不能证明分化关系。feature_genes仅用于拟合后的展示，不参与轨迹构建。
#
# allow_exploratory_low_support = TRUE permits descriptive output when a cell type
# is sparse. Such results are labelled exploratory and are never promoted to a
# replicate-supported Control-to-Case conclusion.
# allow_exploratory_low_support = TRUE允许在细胞不足时生成描述性输出，但结果会被
# 明确标记为探索性，不能升级为具有生物学重复支持的Control到Case结论。
# force_exploratory = TRUE is used when a profile deliberately relaxes internal
# subcluster-size thresholds; it keeps the resulting trajectory exploratory even
# when total cell and replicate counts satisfy the general support thresholds.
# force_exploratory = TRUE用于主动放宽谱系内部亚群细胞数门槛的情况；即使总细胞数
# 和重复数达到通用标准，该轨迹仍会保持探索性标记。
#
# The configured signatures are consensus identity/state panels assembled from
# established liver-cell biology and atlas-level annotations. They are used as
# transparent hypotheses, not as infallible labels. M14B scores them with UCell,
# audits every available gene, and—where an appropriate reference exists—adds an
# independent SingleR reference-label check. Disease-response genes remain in
# feature_genes for visualization unless they are genuinely part of the endpoint
# state definition; this avoids treating inflammation alone as cell identity.
#
# 下列签名是依据经典肝脏细胞生物学及图谱级注释整理的共识身份/状态基因集合，属于
# 可审计的生物学假设，而不是绝对标签。M14B使用UCell进行排名评分，逐一记录实际
# 存在的基因；当存在合适参考时，再加入独立的SingleR参考标签核验。疾病反应基因
# 原则上保留在feature_genes中用于展示，除非它们确实定义该profile的终末状态，
# 从而避免仅凭炎症强弱判断细胞身份。
#
# Method references / 方法依据：
# UCell: https://doi.org/10.1016/j.csbj.2021.06.043
# SingleR: https://doi.org/10.1038/s41590-018-0276-y
# ImmGen: https://doi.org/10.1038/ni.2630
# Liver cell atlas: https://doi.org/10.1038/s41586-021-03973-1
PSEUDOTIME_LINEAGE_PROFILES <- list(
  list(
    priority = 1L,
    profile_id = "hepatocyte_metabolic_state",
    display_name = "Hepatocyte metabolic state",
    target_cell_types = c(
      "Hepatocyte"
    ),
    root_cell_types = c(
      "Hepatocyte"
    ),
    endpoint_cell_types = c(
      "Hepatocyte"
    ),
    feature_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Gys2", "Gbe1", "Ugp2", "Pgm1",
        "Pygl", "Agl", "Gaa", "Gck",
        "Slc2a2", "G6pc", "Pck1", "Irs2"
      )
    } else {
      c(
        "GYS2", "GBE1", "UGP2", "PGM1",
        "PYGL", "AGL", "GAA", "GCK",
        "SLC2A2", "G6PC", "PCK1", "IRS2"
      )
    },
    root_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Gys2", "Gbe1", "Ugp2", "Pgm1",
        "Gck", "Slc2a2", "Irs2"
      )
    } else {
      c(
        "GYS2", "GBE1", "UGP2", "PGM1",
        "GCK", "SLC2A2", "IRS2"
      )
    },
    endpoint_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Pygl", "Agl", "Gaa", "G6pc", "Pck1"
      )
    } else {
      c(
        "PYGL", "AGL", "GAA", "G6PC", "PCK1"
      )
    },
    allow_exploratory_low_support = TRUE,
    min_total_cells = 20L,
    min_cells_per_condition = 2L,
    min_cluster_cells = 5L,
    n_hvg = 1000L,
    n_pcs = 10L,
    cluster_resolution = 0.4
  ),
  list(
    priority = 2L,
    profile_id = "monocyte_macrophage_transition",
    display_name = "Monocyte-Macrophage transition",
    target_cell_types = c(
      "Monocyte",
      "Macrophage"
    ),
    root_cell_types = c(
      "Monocyte"
    ),
    endpoint_cell_types = c(
      "Macrophage"
    ),
    feature_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Ly6c2", "Ccr2", "Lyz2", "Adgre1",
        "Csf1r", "Fcgr1", "Cd68", "Mertk",
        "Apoe", "Trem2", "Spp1", "Il1b"
      )
    } else {
      c(
        "LY6C2", "CCR2", "LYZ", "ADGRE1",
        "CSF1R", "FCGR1A", "CD68", "MERTK",
        "APOE", "TREM2", "SPP1", "IL1B"
      )
    },
    root_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Ly6c2", "Ccr2", "Lyz2"
      )
    } else {
      c(
        "LY6C2", "CCR2", "LYZ"
      )
    },
    endpoint_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Adgre1", "Csf1r", "Fcgr1", "Cd68", "Mertk"
      )
    } else {
      c(
        "ADGRE1", "CSF1R", "FCGR1A", "CD68", "MERTK"
      )
    },
    # ImmGen is a mouse immune-cell reference. These regular-expression
    # patterns are matched against SingleR label.main results. They are used
    # only as auxiliary evidence and never replace the observed annotation.
    # ImmGen是小鼠免疫细胞参考库。以下正则表达式与SingleR的label.main结果
    # 匹配；它们只提供辅助证据，不会覆盖本数据已有的人工/项目注释。
    reference_root_label_patterns = c(
      "^Monocytes?$"
    ),
    reference_endpoint_label_patterns = c(
      "^Macrophages?$"
    ),
    allow_exploratory_low_support = FALSE
  ),
  list(
    priority = 3L,
    profile_id = "kupffer_homeostatic_activated",
    display_name = "Kupffer homeostatic-activated state",
    target_cell_types = c(
      "Kupffer cell"
    ),
    root_cell_types = c(
      "Kupffer cell"
    ),
    endpoint_cell_types = c(
      "Kupffer cell"
    ),
    feature_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Clec4f", "Timd4", "Marco", "Vsig4",
        "Cd5l", "Trem2", "Spp1", "Gpnmb",
        "Lgals3", "Apoe", "Il1b", "Tnf"
      )
    } else {
      c(
        "CLEC4F", "TIMD4", "MARCO", "VSIG4",
        "CD5L", "TREM2", "SPP1", "GPNMB",
        "LGALS3", "APOE", "IL1B", "TNF"
      )
    },
    root_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Clec4f", "Timd4", "Marco", "Vsig4", "Cd5l"
      )
    } else {
      c(
        "CLEC4F", "TIMD4", "MARCO", "VSIG4", "CD5L"
      )
    },
    endpoint_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Trem2", "Spp1", "Gpnmb", "Lgals3",
        "Apoe", "Il1b", "Tnf"
      )
    } else {
      c(
        "TREM2", "SPP1", "GPNMB", "LGALS3",
        "APOE", "IL1B", "TNF"
      )
    },
    allow_exploratory_low_support = TRUE,
    force_exploratory = TRUE,
    min_total_cells = 20L,
    min_cells_per_condition = 2L,
    min_cluster_cells = 5L,
    n_hvg = 1000L,
    n_pcs = 10L,
    cluster_resolution = 0.8
  ),
  list(
    priority = 4L,
    profile_id = "lsec_homeostatic_dysfunctional",
    display_name = "LSEC homeostatic-dysfunctional state",
    target_cell_types = c(
      "LSEC / Endo L"
    ),
    root_cell_types = c(
      "LSEC / Endo L"
    ),
    endpoint_cell_types = c(
      "LSEC / Endo L"
    ),
    feature_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Clec4g", "Stab1", "Stab2", "Lyve1",
        "Kdr", "Klf2", "Plvap", "Cd34",
        "Vwf", "Pecam1", "Col4a1", "Col4a2"
      )
    } else {
      c(
        "CLEC4G", "STAB1", "STAB2", "LYVE1",
        "KDR", "KLF2", "PLVAP", "CD34",
        "VWF", "PECAM1", "COL4A1", "COL4A2"
      )
    },
    root_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Clec4g", "Stab1", "Stab2", "Lyve1",
        "Kdr", "Klf2"
      )
    } else {
      c(
        "CLEC4G", "STAB1", "STAB2", "LYVE1",
        "KDR", "KLF2"
      )
    },
    endpoint_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Plvap", "Cd34", "Vwf", "Pecam1",
        "Col4a1", "Col4a2"
      )
    } else {
      c(
        "PLVAP", "CD34", "VWF", "PECAM1",
        "COL4A1", "COL4A2"
      )
    },
    allow_exploratory_low_support = FALSE,
    min_cluster_cells = 15L
  ),
  list(
    priority = 5L,
    profile_id = "stellate_fibrogenic_activation",
    display_name = "Stellate-fibroblast activation state",
    target_cell_types = c(
      "Fibroblast / stellate cell"
    ),
    root_cell_types = c(
      "Fibroblast / stellate cell"
    ),
    endpoint_cell_types = c(
      "Fibroblast / stellate cell"
    ),
    feature_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Lrat", "Reln", "Des", "Rgs5",
        "Colec10", "Col1a1", "Col1a2", "Col3a1",
        "Timp1", "Acta2", "Tagln", "Postn"
      )
    } else {
      c(
        "LRAT", "RELN", "DES", "RGS5",
        "COLEC10", "COL1A1", "COL1A2", "COL3A1",
        "TIMP1", "ACTA2", "TAGLN", "POSTN"
      )
    },
    root_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Lrat", "Reln", "Des", "Rgs5", "Colec10"
      )
    } else {
      c(
        "LRAT", "RELN", "DES", "RGS5", "COLEC10"
      )
    },
    endpoint_marker_genes = if (
      tolower(SPECIES) == "mouse"
    ) {
      c(
        "Col1a1", "Col1a2", "Col3a1", "Timp1",
        "Acta2", "Tagln", "Postn"
      )
    } else {
      c(
        "COL1A1", "COL1A2", "COL3A1", "TIMP1",
        "ACTA2", "TAGLN", "POSTN"
      )
    },
    allow_exploratory_low_support = FALSE
  )
)

# Default settings are used whenever a profile does not declare an override.
# 某条谱系未单独指定参数时，使用以下通用默认值。
PSEUDOTIME_TARGET_CELL_TYPES <- character(0)
PSEUDOTIME_ROOT_CELL_TYPES <- character(0)
PSEUDOTIME_ENDPOINT_CELL_TYPES <- character(0)
PSEUDOTIME_FEATURE_GENES <- character(0)

# Subclustering parameters used only inside the selected lineage. Recomputing PCA,
# neighbors and UMAP within the subset reveals within-lineage structure that may be
# hidden by global differences among unrelated cell types.
# 仅用于目标谱系内部重新聚类的参数。对子集重新计算PCA、邻居图和UMAP，可以显示
# 在全局不同cell type差异中可能被掩盖的谱系内部结构。
PSEUDOTIME_N_HVG <- 2000L
PSEUDOTIME_N_PCS <- 20L
PSEUDOTIME_CLUSTER_RESOLUTION <- 0.4

# Safety thresholds. The module stops rather than drawing an unstable trajectory
# when the lineage, either condition or eligible subclusters contain too few cells.
# 安全阈值。当目标谱系、任一condition或合格亚群细胞数不足时，模块会停止，而不是
# 强行绘制不稳定轨迹。
PSEUDOTIME_MIN_TOTAL_CELLS <- 100L
PSEUDOTIME_MIN_CELLS_PER_CONDITION <- 20L
PSEUDOTIME_MIN_CLUSTER_CELLS <- 20L

# Minimum fraction of biological replicates from EACH condition that must
# contribute at least one cell to a candidate root or endpoint subcluster. This
# prevents a small, sample-specific population from defining trajectory direction.
# 候选起点或终点亚群在每个condition中必须覆盖的最少生物学重复比例。该阈值可以
# 防止由单一样本或少数样本主导的小亚群决定整条轨迹方向。
PSEUDOTIME_MIN_REPLICATE_COVERAGE <- 0.67

# Evidence weights used to select the root and endpoint subclusters. UCell is
# the primary evidence because its rank-based score is less sensitive than a raw
# mean to library size and a few highly expressed genes. SingleR/ImmGen provides
# an independent reference vote only for compatible mouse immune profiles.
# CONTROL_GROUP/CASE_GROUP remains an explicit directional anchor, but its 25%
# weight prevents the experimental label from overriding contradictory biology.
# If reference evidence is unavailable or inappropriate, available weights are
# renormalized automatically: signature = 2/3 and condition = 1/3. Thus Control/
# Case always influences orientation without being treated as observed time.
#
# 用于选择起点和终点亚群的证据权重。UCell基于基因排名，对文库大小及少数高表达
# 基因的敏感性低于简单表达均值，因此作为主要证据。SingleR/ImmGen只在相容的小鼠
# 免疫谱系中提供独立参考票。CONTROL_GROUP/CASE_GROUP仍是明确的方向锚点，但25%
# 的权重可避免实验分组压过相反的生物学证据。若外部参考不可用或不适用，可用权重
# 会自动重新归一化为signature 2/3、condition 1/3；因此Control/Case始终参与定向，
# 但不会被误解为真实采样时间。
PSEUDOTIME_SELECTION_WEIGHT_SIGNATURE <- 0.50
PSEUDOTIME_SELECTION_WEIGHT_REFERENCE <- 0.25
PSEUDOTIME_SELECTION_WEIGHT_CONDITION <- 0.25

# UCell ranks genes within each cell and evaluates whether a configured signature
# is enriched near the top of that ranking. maxRank limits the evaluated rank
# range; 1500 is the package default and is capped automatically by gene count.
# UCell先在每个细胞内对基因排序，再判断目标签名是否富集于排名前端。maxRank限制
# 参与评价的排名范围；1500为软件包默认值，并会自动受实际基因数上限约束。
PSEUDOTIME_UCELL_MAX_RANK <- 1500L

# Enable the independent ImmGen check when the current species is mouse and the
# selected profile supplies reference label patterns. Reference retrieval or
# matching failure is recorded and causes a documented fallback, not module loss.
# 当物种为mouse且profile提供参考标签模式时，启用独立ImmGen核验。参考数据获取或
# 标签匹配失败会被记录并触发有说明的回退，不会导致整条谱系分析丢失。
PSEUDOTIME_USE_IMMGEN_REFERENCE <- TRUE

# Minimum increase in cell-level median scaled pseudotime required before the
# requested CONTROL_GROUP-to-CASE_GROUP orientation is labelled directionally
# consistent. This is a descriptive effect-size check, not a significance test.
# 将CONTROL_GROUP到CASE_GROUP方向标记为一致所需的cell层面拟时中位数最小增幅。
# 这是描述性效应量检查，不是显著性检验。
PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE <- 0.05

# Number of highly pseudotime-associated genes displayed in the exploratory
# heatmap. This heatmap is descriptive and is not a branch-specific significance
# test equivalent to Monocle BEAM.
# 探索性拟时热图展示的高关联基因数量。该热图属于描述性结果，不能等同于Monocle
# BEAM的分支特异性显著性检验。
PSEUDOTIME_HEATMAP_TOP_GENES <- 30L
PSEUDOTIME_HEATMAP_BINS <- 40L

# Reproducibility / 可重复性
# Fixed seed used by stochastic algorithms such as UMAP. Reproducibility also
# depends on software versions, numerical libraries and parallel configuration,
# all of which should be retained with the environment log.
# 随机算法（如UMAP）使用的固定种子。完整可重复性还取决于软件版本、数值计算库和
# 并行设置，因此应同时保存环境日志。
RANDOM_SEED <- 20260811
set.seed(RANDOM_SEED)


# ============================================================
# 0A. SPECIES CONFIGURATION / 物种配置
# ============================================================

# Centralized species-specific regular expressions used to identify mitochondrial
# and ribosomal genes from standard gene symbols. These patterns are validated in
# preflight checks; custom annotations may require explicit adjustment.
# 集中定义依据标准基因符号识别线粒体和核糖体基因的物种特异性正则表达式。预检会
# 验证这些模式；若使用自定义基因注释，可能需要明确调整。

SPECIES_CONFIG <- list(
  mouse = list(
    mito_pattern = "^mt-",
    ribo_pattern = "^Rp[sl]",
    orgdb_package = "org.Mm.eg.db",
    kegg_code = "mmu",
    cellchat_database = "mouse"
  ),
  human = list(
    mito_pattern = "^MT-",
    ribo_pattern = "^RP[SL]",
    orgdb_package = "org.Hs.eg.db",
    kegg_code = "hsa",
    cellchat_database = "human"
  ),
  rat = list(
    mito_pattern = "^mt-",
    ribo_pattern = "^Rp[sl]",
    orgdb_package = "org.Rn.eg.db",
    kegg_code = "rno",
    cellchat_database = NA_character_
  )
)

species_key <- tolower(SPECIES)

if (!species_key %in% names(SPECIES_CONFIG)) {
  stop(
    "Unsupported SPECIES / 暂不支持该物种: ", SPECIES,
    "\nSupported / 当前支持: ",
    paste(names(SPECIES_CONFIG), collapse = ", "),
    call. = FALSE
  )
}

species_cfg <- SPECIES_CONFIG[[species_key]]


# ============================================================
# 0B. INTERNAL MODULE SWITCHES / 内部模块开关
# ============================================================

# These switches are derived automatically from RUN_FROM, RUN_TO and the optional
# branch settings. They document execution state and should not normally be edited
# by hand; changing the user-facing module range above is safer and auditable.
# 以下开关由RUN_FROM、RUN_TO和可选分支设置自动生成，用于记录本次执行状态，通常
# 不应手动修改；调整上方用户配置中的模块范围更安全，也更便于审计。

run_step <- function(x) RUN_FROM <= x && RUN_TO >= x

RUN_M00_ENVIRONMENT        <- run_step(0)
RUN_M01_SAMPLE_METADATA    <- run_step(1)
RUN_M02_ACQUIRE_DATA       <- run_step(2)
RUN_M03_INPUT_AUDIT        <- run_step(3)
RUN_M04_READ_INPUT         <- run_step(4)
RUN_M05_DEMULTIPLEX        <- run_step(5)
RUN_M06_TISSUE_SUBSET      <- run_step(6)
RUN_M07_QC                 <- run_step(7)
RUN_M08_DOUBLETS           <- run_step(8)
RUN_M09_NORMALIZE_HVG      <- run_step(9)
RUN_M10_PCA_CLUSTER_UMAP   <- run_step(10)
RUN_M10B_INTEGRATION       <- RUN_OPTIONAL_INTEGRATION && run_step(10)
RUN_M11_MARKERS_ANNOTATION <- run_step(11)
RUN_M12_COMPOSITION        <- run_step(12)
RUN_M13_PSEUDOBULK_DE      <- RUN_PSEUDOBULK && run_step(13)
RUN_M13B_FUNCTIONAL        <- RUN_PSEUDOBULK &&
  RUN_FUNCTIONAL_INTERPRETATION && run_step(13)
RUN_M13C_GSVA              <- RUN_PSEUDOBULK &&
  RUN_GSVA && run_step(13)
RUN_M14_FOCUSED_ANALYSIS   <- run_step(14)
RUN_M14B_PSEUDOTIME        <- RUN_PSEUDOTIME && run_step(14)
RUN_M14C_CELLCHAT          <- RUN_CELLCHAT && run_step(14)
RUN_M15_FINAL_SAVE         <- run_step(15)

DATA_SOURCE <- if (USE_LOCAL_FILES) "LOCAL" else "GEO"
GSE_ID <- if (USE_LOCAL_FILES) "" else DATASET_ID
DOWNLOAD_GEO_DATA <- !USE_LOCAL_FILES
FORCE_REEXTRACT <- FALSE
HAS_HTO <- identical(MULTIPLEXING_MODE, "HTO")

# ============================================================
# 0C. PROJECT PATHS / 项目路径
# ============================================================

# All result paths are anchored to the current R project directory. The workflow
# separates configuration, raw input, processed checkpoints and module-specific
# results so that source data and derived outputs remain distinguishable.
# 所有结果路径均以当前R项目目录为基准。流程将配置、原始输入、处理后检查点和
# 模块结果分别存放，确保源数据与派生结果能够清楚区分。

PROJECT_DIR <- normalizePath(
  ".",
  winslash = "/",
  mustWork = TRUE
)

DIR_CONFIG <- file.path(PROJECT_DIR, "config")
DIR_RAW <- file.path(PROJECT_DIR, "data", "raw")
DIR_PROCESSED <- file.path(PROJECT_DIR, "data", "processed", DATASET_ID)
DIR_RESULTS <- file.path(PROJECT_DIR, "results", DATASET_ID)
DIR_MODULES <- file.path(DIR_RESULTS, "modules")

# INPUT_DIR may be relative to the RStudio Project or an absolute path.
# INPUT_DIR可以是相对于RStudio Project的路径，也可以是绝对路径。
is_absolute_path <- grepl("^(/|[A-Za-z]:[/\\\\])", INPUT_DIR)

LOCAL_INPUT_DIR <- if (is_absolute_path) {
  normalizePath(INPUT_DIR, winslash = "/", mustWork = FALSE)
} else {
  file.path(PROJECT_DIR, INPUT_DIR)
}

all_base_dirs <- c(
  DIR_CONFIG,
  DIR_RAW,
  DIR_PROCESSED,
  DIR_RESULTS,
  DIR_MODULES
)

invisible(lapply(
  all_base_dirs,
  dir.create,
  recursive = TRUE,
  showWarnings = FALSE
))

if (!nzchar(SAMPLE_METADATA_FILE)) {
  SAMPLE_METADATA_FILE <- file.path(
    DIR_CONFIG,
    paste0(DATASET_ID, "_sample_metadata.csv")
  )
}

cat(
  "\nDataset / 数据集: ", DATASET_ID,
  "\nInput directory / 输入目录: ", LOCAL_INPUT_DIR,
  "\nModules this run / 本次模块: M", sprintf("%02d", RUN_FROM),
  " to M", sprintf("%02d", RUN_TO), "\n\n",
  sep = ""
)


# ============================================================
# 0D. MODULE HELPERS / 模块辅助函数
# ============================================================

# Helper functions below provide consistent directory creation, checkpoint loading,
# metadata validation, plotting and logging across modules. They are implementation
# infrastructure rather than dataset-specific settings and should be changed only
# when the workflow itself is deliberately revised and revalidated.
# 以下辅助函数用于在各模块之间统一创建目录、读取检查点、验证metadata、绘图和记录
# 日志。它们属于流程基础设施而非数据集配置，只有在明确修改并重新验证整个工作流
# 时才应调整。

module_dir <- function(module_name) {

  x <- file.path(
    DIR_MODULES,
    module_name
  )

  dirs <- c(
    x,
    file.path(x, "tables"),
    file.path(x, "figures"),
    file.path(x, "logs")
  )

  invisible(lapply(
    dirs,
    dir.create,
    recursive = TRUE,
    showWarnings = FALSE
  ))

  x
}


checkpoint_file <- function(name) {
  file.path(
    DIR_PROCESSED,
    paste0(name, ".rds")
  )
}


save_checkpoint <- function(object, name) {

  path <- checkpoint_file(name)

  saveRDS(
    object,
    path
  )

  message(
    "Checkpoint saved: ",
    path
  )

  invisible(path)
}


# ============================================================
# TRANSPARENT FIGURE OUTPUT / 透明图片输出
# ============================================================
#
# PNG device transparency alone does not remove the backgrounds drawn by
# ggplot2 themes. The helper below therefore makes the plot canvas, plotting
# panels, legends, legend keys and facet strips transparent before saving.
# It supports both ordinary ggplot objects and Seurat/patchwork composite plots.
#
# 仅设置PNG图形设备的透明背景，并不能去除ggplot2主题自身绘制的背景。因此，
# 下方辅助函数会在保存前统一将整张图片画布、绘图区、图例、图例键和分面标签背景
# 设为透明。该函数同时支持普通ggplot对象和Seurat/patchwork组合图。
#
# This changes figure presentation only. It does not modify expression values,
# dimensional reductions, clustering, annotation or statistical analysis.
#
# 该设置仅改变图片呈现方式，不会修改表达量、降维、聚类、注释或统计分析结果。
# ============================================================

make_plot_background_transparent <- function(plot) {

  transparent_theme <- ggplot2::theme(
    plot.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    panel.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    legend.key = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    ),
    strip.background = ggplot2::element_rect(
      fill = "transparent",
      colour = NA
    )
  )

  if (inherits(plot, "patchwork")) {
    plot & transparent_theme
  } else {
    plot + transparent_theme
  }
}


save_plot_transparent <- function(
    filename,
    plot,
    ...
) {

  plot <- make_plot_background_transparent(
    plot
  )

  ggplot2::ggsave(
    filename = filename,
    plot = plot,
    ...,
    bg = "transparent"
  )
}


load_checkpoint <- function(name) {

  path <- checkpoint_file(name)

  if (!file.exists(path)) {
    stop(
      "Checkpoint不存在：\n",
      path,
      "\n请先完成上一模块。"
    )
  }

  readRDS(path)
}


# ============================================================
# SEURAT v5 RNA-LAYER COMPATIBILITY
# Seurat v5 RNA layer兼容处理
# ============================================================
#
# Seurat v5 Assay5 can contain multiple layers after objects
# from different libraries have been merged, e.g. counts.1,
# counts.2, counts.3. Some Bioconductor conversions and legacy
# data-access functions expect one coherent count layer.
#
# Seurat v5的Assay5在多个library对象合并后可能保留多个layers，
# 例如counts.1、counts.2、counts.3。部分Bioconductor转换函数
# 和传统数据读取函数要求一个统一的count layer。
#
# JoinLayers() is NOT batch correction, integration, normalization,
# or biological pooling. It only joins split assay layers so that
# downstream functions can access the matrix consistently.
#
# JoinLayers()不是批次校正、integration、normalization，
# 也不会把生物学重复混为一个样本。sample_id、
# biological_replicate等cell metadata都会保留。
# ============================================================

join_assay_layers_if_needed <- function(
    object,
    assay = "RNA",
    verbose = TRUE
) {

  if (!(assay %in% SeuratObject::Assays(object))) {
    return(object)
  }

  assay_object <- object[[assay]]

  if (!inherits(assay_object, "Assay5")) {
    return(object)
  }

  layer_names <- SeuratObject::Layers(
    assay_object
  )

  split_layer_pattern <- "^(counts|data|scale\\.data)\\."
  has_split_layers <- any(
    grepl(
      split_layer_pattern,
      layer_names
    )
  )

  if (!has_split_layers) {
    return(object)
  }

  if (verbose) {
    message(
      "Seurat v5 multi-layer assay detected / 检测到Seurat v5多layer assay: ",
      assay,
      " [",
      paste(layer_names, collapse = ", "),
      "]"
    )

    message(
      "Joining layers for downstream compatibility / ",
      "正在合并layers以保证后续兼容性。"
    )
  }

  object <- SeuratObject::JoinLayers(
    object,
    assay = assay
  )

  object
}


prepare_for_sce <- function(
    object,
    assay = "RNA"
) {

  object <- join_assay_layers_if_needed(
    object,
    assay = assay,
    verbose = TRUE
  )

  DefaultAssay(object) <- assay

  object
}


write_module_status <- function(
    module_name,
    status,
    message_text = ""
) {

  d <- module_dir(module_name)

  row <- data.frame(
    module = module_name,
    status = status,
    time = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    ),
    message = message_text,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    row,
    file.path(
      d,
      "logs",
      "module_status.csv"
    )
  )
}


safe_stop_module <- function(
    module_name,
    message_text
) {

  write_module_status(
    module_name,
    "STOPPED",
    message_text
  )

  stop(
    message_text,
    call. = FALSE
  )
}


validate_complete_annotation <- function(
    object,
    module_name
) {

  metadata <- object[[]]

  if (!"cell_type" %in% colnames(metadata)) {
    safe_stop_module(
      module_name,
      "cell_type metadata is missing / 缺少cell_type metadata。"
    )
  }

  cell_type_values <- trimws(
    as.character(
      metadata$cell_type
    )
  )

  incomplete <- (
    is.na(cell_type_values) |
      !nzchar(cell_type_values) |
      grepl(
        "^Cluster_",
        cell_type_values
      )
  )

  if (any(incomplete)) {
    safe_stop_module(
      module_name,
      paste0(
        sum(incomplete),
        " cells have missing or placeholder annotations / ",
        "个cell仍为缺失或占位注释。"
      )
    )
  }

  invisible(TRUE)
}


first_existing <- function(paths) {

  paths <- paths[
    !is.na(paths) &
      nzchar(paths)
  ]

  existing <- paths[
    file.exists(paths)
  ]

  if (length(existing) == 0L) {
    return(NA_character_)
  }

  existing[[1L]]
}


clean_file_stem <- function(x) {

  x <- basename(x)

  x <- sub(
    "\\.gz$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- sub(
    "matrix\\.mtx$",
    "",
    x,
    ignore.case = TRUE
  )

  x <- sub(
    "[_.-]+$",
    "",
    x
  )

  x
}


extract_gsm <- function(x) {

  m <- regexpr(
    "GSM[0-9]+",
    x,
    ignore.case = TRUE
  )

  if (m[[1]] == -1L) {
    return(NA_character_)
  }

  toupper(
    regmatches(
      x,
      m
    )
  )
}




# ============================================================
# BEGINNER CORE DICTIONARY / 零基础核心词典
# ============================================================
#
# count / counts:
#   原始表达计数，10x数据中通常是UMI counts。
#   Raw expression counts; usually UMI counts in 10x data.
#
# UMI (Unique Molecular Identifier):
#   用于降低PCR重复扩增造成的计数偏差。
#   Molecular tags used to reduce PCR-duplication bias.
#
# barcode:
#   一个液滴/细胞的标识符。
#   Identifier for a droplet/cell.
#
# sparse matrix:
#   大多数gene × cell位置为0，只存非零值以节省内存。
#   Stores only non-zero entries to save memory.
#
# assay:
#   Seurat对象中的一种数据模态，例如RNA或HTO。
#   A data modality inside Seurat, e.g. RNA or HTO.
#
# metadata:
#   与每个cell关联的信息，如sample、condition、QC、cluster、cell type。
#   Cell-associated information such as sample, condition, QC and annotation.
#
# normalization:
#   降低细胞间测序深度/捕获效率差异的影响。
#   Adjusts for cell-specific sequencing/capture depth.
#
# HVG (Highly Variable Gene):
#   在技术/基础变化之外表现出较强生物学异质性的基因。
#   A gene with substantial biological variability beyond baseline noise.
#
# PCA:
#   将数千基因压缩为较少的主成分。
#   Compresses thousands of genes into a smaller number of components.
#
# KNN/SNN:
#   根据PCA空间中的相似性建立细胞邻居网络。
#   Cell-neighbor graphs based on similarity in PCA space.
#
# cluster:
#   算法定义的表达相似细胞群；编号本身没有生物学含义。
#   Algorithm-defined group of transcriptionally similar cells.
#
# marker:
#   帮助解释某个cluster身份的差异表达基因。
#   Gene that helps distinguish and interpret a cluster.
#
# biological replicate:
#   独立动物/患者/实验单位，是正式组间推断的统计重复。
#   Independent animal/patient/experimental unit for inference.
#
# pseudobulk:
#   同一生物学重复、同一细胞类型内，将raw counts按基因求和。
#   Sum raw counts within one replicate and one cell type.
#
# pseudotime:
#   根据同一谱系内细胞的转录连续性推断的相对状态顺序，不是真实采样时间，也不能
#   单独证明细胞来源、分化方向或处理效应的因果关系。
#   Relative state ordering inferred from transcriptional continuity within one
#   lineage; it is not observed time and does not by itself prove lineage,
#   direction or a causal treatment effect.
#
# logFC:
#   CASE_GROUP相对CONTROL_GROUP的log2 fold change。
#   log2 fold change of CASE_GROUP relative to CONTROL_GROUP.
#
# FDR:
#   多重检验校正后的错误发现率。
#   False-discovery rate after multiple-testing correction.
#
# ============================================================

# ============================================================
# M00. SOFTWARE ENVIRONMENT / 软件环境
#
# 【Module scope / 模块范围】
# This module validates computational prerequisites only; it does not read or
# transform biological data. Its records define the software context in which
# all later numerical results were generated.
# 本模块仅验证计算环境，不读取或转换生物学数据。其输出用于说明后续所有数值结果
# 所依赖的软件环境，是复现和审阅分析的起点。
#
# 【Decision boundary / 决策边界】
# Package installation changes the local software environment but not the analysis
# profile. If exact version reproduction is required, use a locked environment and
# set INSTALL_PACKAGES = FALSE after dependencies have been prepared.
# 软件包安装会改变本地软件环境，但不改变分析配置。若需要严格复现软件版本，应先
# 建立锁定环境，再将INSTALL_PACKAGES设为FALSE。
#
# 【Purpose / 目的】
# 检查并安装所需R包，记录版本、随机种子和关键分析参数。
# Check/install required packages and record package versions,
# random seed and key analysis settings.
#
# 【Input / 输入】
# 第0节用户参数；不读取表达矩阵。
# User settings from Section 0; no expression matrix is read here.
#
# 【Output / 输出】
#   tables/package_versions.csv
#   tables/analysis_parameters.csv
#   logs/sessionInfo.txt
#
# 【Why it matters / 为什么重要】
# 单细胞分析结果可能随Seurat/Bioconductor版本变化。
# Single-cell results may change across package versions.
# 记录环境可以让未来自己或他人复现同一分析。
# Recording the environment supports reproducibility.
#
# 【Check after running / 运行后检查】
# 1. 所有包都有版本号 / every package has a version.
# 2. 没有missing package / no package is missing.
# 3. DATASET_ID、CONTROL_GROUP、CASE_GROUP等设置正确。
#    Confirm dataset and group settings.
# ============================================================

if (RUN_M00_ENVIRONMENT) {

  MODULE <- "M00_environment"
  MDIR <- module_dir(MODULE)

  options(timeout = 3600)
  options(
    future.globals.maxSize =
      20 * 1024^3
  )

  cran_packages <- c(
    "Seurat",
    "SeuratObject",
    "ggplot2",
    "dplyr",
    "tidyr",
    "tibble",
    "data.table",
    "patchwork",
    "ggrepel",
    "Matrix",
    "future",
    "scales",
    "clustree",
    "ragg",
    "circlize",
    "renv"
  )

  if (DATA_SOURCE == "GEO") {
    cran_packages <- unique(
      c(
        cran_packages,
        "rentrez"
      )
    )
  }

  if (RUN_PSEUDOTIME) {
    cran_packages <- unique(
      c(
        cran_packages,
        "mgcv"
      )
    )
  }

  if (RUN_CELLCHAT) {
    cran_packages <- unique(
      c(
        cran_packages,
        "remotes",
        "igraph",
        "ggraph",
        "tidygraph"
      )
    )
  }

  bioc_packages <- c(
    "SingleCellExperiment",
    "SummarizedExperiment",
    "S4Vectors",
    "scuttle",
    "scrapper",
    "scran",
    "scDblFinder",
    "edgeR",
    "limma",
    "BiocParallel",
    "ComplexHeatmap"
  )

  if (RUN_PSEUDOTIME) {
    bioc_packages <- unique(
      c(
        bioc_packages,
        "slingshot",
        "DelayedMatrixStats",
        "UCell"
      )
    )

    if (
      isTRUE(
        PSEUDOTIME_USE_IMMGEN_REFERENCE
      ) &&
      identical(
        species_key,
        "mouse"
      )
    ) {
      bioc_packages <- unique(
        c(
          bioc_packages,
          "SingleR",
          "celldex"
        )
      )
    }
  }

  if (RUN_FUNCTIONAL_INTERPRETATION) {
    bioc_packages <- unique(
      c(
        bioc_packages,
        "AnnotationDbi",
        "clusterProfiler",
        "enrichplot",
        "ggkegg",
        species_cfg$orgdb_package
      )
    )
  }

  if (RUN_GSVA) {
    bioc_packages <- unique(
      c(
        bioc_packages,
        "AnnotationDbi",
        "clusterProfiler",
        "GSVA",
        "KEGGREST",
        species_cfg$orgdb_package
      )
    )
  }

  if (DATA_SOURCE == "GEO") {
    bioc_packages <- unique(
      c(
        bioc_packages,
        "GEOquery"
      )
    )
  }

  install_missing_cran <- function(packages) {

    missing <- packages[
      !vapply(
        packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
      )
    ]

    if (length(missing) > 0L) {
      install.packages(
        missing,
        dependencies = TRUE
      )
    }
  }


  install_missing_bioc <- function(packages) {

    if (
      !requireNamespace(
        "BiocManager",
        quietly = TRUE
      )
    ) {
      install.packages(
        "BiocManager"
      )
    }

    missing <- packages[
      !vapply(
        packages,
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
      )
    ]

    if (length(missing) > 0L) {
      BiocManager::install(
        missing,
        ask = FALSE,
        update = FALSE
      )
    }
  }


  if (INSTALL_PACKAGES) {
    install_missing_cran(
      cran_packages
    )
    install_missing_bioc(
      bioc_packages
    )

    # CellChat is maintained in its official GitHub repository rather than the
    # current CRAN/Bioconductor release channels. Installation is therefore
    # explicit and occurs only when the optional communication branch is enabled.
    # CellChat当前由官方GitHub仓库维护，而非通过当前CRAN/Bioconductor
    # 发行渠道提供。因此只在启用可选细胞互作分支时显式安装。
    if (
      RUN_CELLCHAT &&
      !requireNamespace(
        "CellChat",
        quietly = TRUE
      )
    ) {
      remotes::install_github(
        "jinworks/CellChat",
        upgrade = "never",
        dependencies = TRUE
      )
    }
  }

  all_packages <- unique(
    c(
      cran_packages,
      bioc_packages,
      if (RUN_CELLCHAT) "CellChat" else character(0)
    )
  )

  available <- vapply(
    all_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )

  if (any(!available)) {
    safe_stop_module(
      MODULE,
      paste0(
        "以下软件包未安装成功：",
        paste(
          all_packages[!available],
          collapse = ", "
        )
      )
    )
  }

  if (USE_PARALLEL) {
    future::plan(
      future::multisession,
      workers = N_WORKERS
    )
  } else {
    future::plan("sequential")
  }

  versions <- data.frame(
    package = all_packages,
    version = vapply(
      all_packages,
      function(pkg) {
        as.character(
          utils::packageVersion(pkg)
        )
      },
      FUN.VALUE = character(1)
    ),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    versions,
    file.path(
      MDIR,
      "tables",
      "package_versions.csv"
    )
  )

  writeLines(
    capture.output(
      sessionInfo()
    ),
    file.path(
      MDIR,
      "logs",
      "sessionInfo.txt"
    )
  )

  pseudotime_profile_ids <- vapply(
    PSEUDOTIME_LINEAGE_PROFILES,
    function(profile) {
      as.character(
        profile$profile_id
      )
    },
    character(1)
  )

  pseudotime_profile_target_cell_types <- vapply(
    PSEUDOTIME_LINEAGE_PROFILES,
    function(profile) {
      paste(
        profile$target_cell_types,
        collapse = "|"
      )
    },
    character(1)
  )

  parameters <- data.frame(
    parameter = c(
      "DATASET_ID",
      "DATA_SOURCE",
      "INPUT_FORMAT",
      "GSE_ID",
      "CONTROL_GROUP",
      "CASE_GROUP",
      "TARGET_TISSUE",
      "TISSUE_ASSIGNMENT_MODE",
      "MULTIPLEXING_MODE",
      "QC_MODE",
      "NORMALIZATION_METHOD",
      "N_PCS_USE",
      "CLUSTER_RESOLUTION",
      "MIN_CELLS_PER_PSEUDOBULK",
      "PSEUDOBULK_REPLICATE_COL",
      "RUN_FUNCTIONAL_INTERPRETATION",
      "DEG_CLASS_HEATMAP_TOP_PER_CLASS",
      "ENRICHMENT_MIN_MAPPED_DEGS",
      "ENRICHMENT_FDR_THRESHOLD",
      "ENRICHMENT_TOP_TERMS_PER_CLASS",
      "RUN_GGKEGG",
      "GGKEGG_PATHWAYS_PER_CLASS",
      "GGKEGG_MAX_PATHWAYS",
      "GGKEGG_MAX_GENES_IN_COMPANION",
      "RUN_GSVA",
      "GSVA_MIN_GENE_SET_SIZE",
      "GSVA_MAX_GENE_SET_SIZE",
      "GSVA_LOGCPM_PRIOR_COUNT",
      "GSVA_FDR_THRESHOLD",
      "GSVA_MIN_ABS_SCORE_DIFFERENCE",
      "GSVA_TOP_PATHWAYS_PER_CELL_TYPE",
      "GSVA_OVERVIEW_MAX_PATHWAYS",
      "RUN_CELLCHAT",
      "CELLCHAT_MIN_CELLS_PER_TYPE_CONDITION",
      "CELLCHAT_MIN_REPLICATES_PER_TYPE_CONDITION",
      "CELLCHAT_MIN_CELLS_FOR_INTERACTION",
      "CELLCHAT_TOP_LR_PAIRS",
      "CELLCHAT_NETWORK_TOP_EDGES",
      "LINEAGE_SUBCLUSTER_TOP_MARKERS",
      "RUN_PSEUDOTIME",
      "PSEUDOTIME_PROFILE_COUNT",
      "PSEUDOTIME_PROFILE_IDS",
      "PSEUDOTIME_PROFILE_TARGET_CELL_TYPES",
      "PSEUDOTIME_SELECTION_WEIGHT_SIGNATURE",
      "PSEUDOTIME_SELECTION_WEIGHT_REFERENCE",
      "PSEUDOTIME_SELECTION_WEIGHT_CONDITION",
      "PSEUDOTIME_UCELL_MAX_RANK",
      "PSEUDOTIME_USE_IMMGEN_REFERENCE",
      "RANDOM_SEED"
    ),
    value = as.character(c(
      DATASET_ID,
      DATA_SOURCE,
      INPUT_FORMAT,
      GSE_ID,
      CONTROL_GROUP,
      CASE_GROUP,
      TARGET_TISSUE,
      TISSUE_ASSIGNMENT_MODE,
      MULTIPLEXING_MODE,
      QC_MODE,
      NORMALIZATION_METHOD,
      N_PCS_USE,
      CLUSTER_RESOLUTION,
      MIN_CELLS_PER_PSEUDOBULK,
      PSEUDOBULK_REPLICATE_COL,
      RUN_FUNCTIONAL_INTERPRETATION,
      DEG_CLASS_HEATMAP_TOP_PER_CLASS,
      ENRICHMENT_MIN_MAPPED_DEGS,
      ENRICHMENT_FDR_THRESHOLD,
      ENRICHMENT_TOP_TERMS_PER_CLASS,
      RUN_GGKEGG,
      GGKEGG_PATHWAYS_PER_CLASS,
      GGKEGG_MAX_PATHWAYS,
      GGKEGG_MAX_GENES_IN_COMPANION,
      RUN_GSVA,
      GSVA_MIN_GENE_SET_SIZE,
      GSVA_MAX_GENE_SET_SIZE,
      GSVA_LOGCPM_PRIOR_COUNT,
      GSVA_FDR_THRESHOLD,
      GSVA_MIN_ABS_SCORE_DIFFERENCE,
      GSVA_TOP_PATHWAYS_PER_CELL_TYPE,
      GSVA_OVERVIEW_MAX_PATHWAYS,
      RUN_CELLCHAT,
      CELLCHAT_MIN_CELLS_PER_TYPE_CONDITION,
      CELLCHAT_MIN_REPLICATES_PER_TYPE_CONDITION,
      CELLCHAT_MIN_CELLS_FOR_INTERACTION,
      CELLCHAT_TOP_LR_PAIRS,
      CELLCHAT_NETWORK_TOP_EDGES,
      LINEAGE_SUBCLUSTER_TOP_MARKERS,
      RUN_PSEUDOTIME,
      length(
        PSEUDOTIME_LINEAGE_PROFILES
      ),
      paste(
        pseudotime_profile_ids,
        collapse = ";"
      ),
      paste(
        pseudotime_profile_target_cell_types,
        collapse = ";"
      ),
      PSEUDOTIME_SELECTION_WEIGHT_SIGNATURE,
      PSEUDOTIME_SELECTION_WEIGHT_REFERENCE,
      PSEUDOTIME_SELECTION_WEIGHT_CONDITION,
      PSEUDOTIME_UCELL_MAX_RANK,
      PSEUDOTIME_USE_IMMGEN_REFERENCE,
      RANDOM_SEED
    )),
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    parameters,
    file.path(
      MDIR,
      "tables",
      "analysis_parameters.csv"
    )
  )

  write_module_status(
    MODULE,
    "OK"
  )
}

# ============================================================
# RUNTIME PACKAGE ATTACHMENT / 运行时软件包加载
# ============================================================
# This block is outside M00 so checkpoint restarts work in a clean R session.
# 本代码块位于M00之外，确保从任意checkpoint启动时都能加载运行时依赖。
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(data.table)
  library(patchwork)
  library(Matrix)
  library(clustree)
})

# ============================================================
# M01. SAMPLE METADATA / 样本信息
#
# 【Module scope / 模块范围】
# This module establishes the experimental-design table that links technical input
# units to independent biological specimens and conditions. No expression-based
# inference is used to repair or guess missing sample identities.
# 本模块建立实验设计表，将技术输入单位对应到独立生物学样本和实验分组。流程不会
# 根据表达数据推测或自动补全缺失的样本身份。
#
# 【Quality gate / 质量门槛】
# Every sample key must be unique, every replicate must map unambiguously to one
# condition, and both comparison groups must contain genuinely independent
# biological replicates. Technical libraries from the same specimen remain the
# same replicate unless the experimental design states otherwise.
# 每个样本键必须唯一，每个生物学重复必须明确对应一个condition，且两组均需包含
# 真正独立的生物学重复。同一标本产生的技术文库仍属于同一重复，除非实验设计另有
# 明确规定。
#
# 【Purpose / 目的】
# 建立“文件/文库 → 生物学重复 → 实验组 → 组织”的对应关系。
# Map each input library to its biological replicate, condition and tissue.
#
# 【Input / 输入】
# config/<DATASET_ID>_sample_metadata.csv
#
# 【最少需要的列 / Minimum required columns】
#   sample_id
#   biological_replicate
#   condition
#   tissue
#
# 【Why it matters / 为什么重要】
# 单细胞统计最常见错误之一，是把细胞数或技术library数当成生物学重复。
# A common single-cell mistake is treating cells or technical libraries
# as independent biological replicates.
#
# pseudobulk最终按biological_replicate而不是cell运行。
# Pseudobulk inference will use biological_replicate, not individual cells.
#
# 【Output / 输出】
#   sample_metadata_used.csv
#   biological_replicates_by_condition.csv
#   M01_sample_metadata.rds
#
# 【Check after running / 运行后检查】
# Control和Case各有多少只独立动物/患者？
# How many independent animals/patients are present in each condition?
# ============================================================

if (RUN_M01_SAMPLE_METADATA) {

  MODULE <- "M01_sample_metadata"
  MDIR <- module_dir(MODULE)

  if (!file.exists(SAMPLE_METADATA_FILE)) {

    template <- data.frame(
      sample_id = c(
        "Sample_1",
        "Sample_2",
        "Sample_3",
        "Sample_4"
      ),
      biological_replicate = c(
        "Animal_1",
        "Animal_2",
        "Animal_3",
        "Animal_4"
      ),
      condition = rep(
        c(
          CONTROL_GROUP,
          CASE_GROUP
        ),
        each = 2L
      ),
      tissue = rep(
        TARGET_TISSUE,
        4
      ),
      sex = NA_character_,
      batch = NA_character_,
      input_path = NA_character_,
      stringsAsFactors = FALSE
    )

    colnames(template)[1:4] <- c(
      SAMPLE_ID_COL,
      BIOLOGICAL_REPLICATE_COL,
      CONDITION_COL,
      TISSUE_COL
    )

    data.table::fwrite(
      template,
      SAMPLE_METADATA_FILE
    )

    safe_stop_module(
      MODULE,
      paste0(
        "已生成sample metadata模板：\n",
        SAMPLE_METADATA_FILE,
        "\n请先按真实实验设计填写，再重新运行M01。"
      )
    )
  }


  sample_metadata <- data.table::fread(
    SAMPLE_METADATA_FILE,
    data.table = FALSE
  )

  required_cols <- c(
    SAMPLE_ID_COL,
    BIOLOGICAL_REPLICATE_COL,
    CONDITION_COL,
    TISSUE_COL
  )

  missing_cols <- setdiff(
    required_cols,
    colnames(sample_metadata)
  )

  if (length(missing_cols) > 0L) {
    safe_stop_module(
      MODULE,
      paste0(
        "sample_metadata缺少必要列：",
        paste(
          missing_cols,
          collapse = ", "
        )
      )
    )
  }


  sample_ids <- trimws(
    as.character(
      sample_metadata[[SAMPLE_ID_COL]]
    )
  )

  if (any(is.na(sample_ids) | !nzchar(sample_ids))) {
    safe_stop_module(
      MODULE,
      "sample_id contains missing/empty values / sample_id存在缺失或空值。"
    )
  }

  if (anyDuplicated(sample_ids) > 0L) {
    safe_stop_module(
      MODULE,
      paste0(
        "Duplicated sample IDs / 重复sample ID: ",
        paste(
          unique(
            sample_ids[duplicated(sample_ids)]
          ),
          collapse = ", "
        )
      )
    )
  }

  replicate_condition_n <- sample_metadata |>
    dplyr::distinct(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]]
    ) |>
    dplyr::count(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      name = "n_conditions"
    )

  if (any(replicate_condition_n$n_conditions != 1L)) {
    safe_stop_module(
      MODULE,
      paste0(
        "A biological replicate maps to multiple conditions / ",
        "同一biological replicate对应多个condition。"
      )
    )
  }


  # 不允许缺失真正的生物学重复。
  if (
    any(
      is.na(
        sample_metadata[[BIOLOGICAL_REPLICATE_COL]]
      )
    ) ||
    any(
      trimws(
        as.character(
          sample_metadata[[BIOLOGICAL_REPLICATE_COL]]
        )
      ) == ""
    )
  ) {
    safe_stop_module(
      MODULE,
      "biological_replicate存在空值。必须明确每个文库来自哪只独立动物/患者。"
    )
  }


  # 检查每个condition的独立重复数。
  replicate_summary <- sample_metadata |>
    dplyr::distinct(
      .data[[CONDITION_COL]],
      .data[[BIOLOGICAL_REPLICATE_COL]]
    ) |>
    dplyr::count(
      .data[[CONDITION_COL]],
      name = "n_biological_replicates"
    )

  data.table::fwrite(
    sample_metadata,
    file.path(
      MDIR,
      "tables",
      "sample_metadata_used.csv"
    )
  )

  data.table::fwrite(
    replicate_summary,
    file.path(
      MDIR,
      "tables",
      "biological_replicates_by_condition.csv"
    )
  )

  print(sample_metadata)
  print(replicate_summary)

  save_checkpoint(
    sample_metadata,
    "M01_sample_metadata"
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M02. DATA ACQUISITION / 获取数据
#
# 【Module scope / 模块范围】
# This module locates or downloads source files and creates a file manifest. It
# does not assume that a file labelled "processed" contains normalized values;
# matrix contents are evaluated in later modules.
# 本模块仅定位或下载源文件并生成文件清单，不会因为文件标记为"processed"就假定
# 其中是标准化表达值；矩阵内容将在后续模块中进一步检查。
#
# 【Provenance requirement / 来源记录要求】
# Retain accession, download location and original filenames. For local data,
# INPUT_DIR should point to an immutable source or a documented prepared object so
# that the analytical input can be reconstructed independently of result folders.
# 应保留登录号、下载位置和原始文件名。本地数据的INPUT_DIR应指向不可变的源文件
# 或具有说明的准备对象，以便在不依赖结果目录的情况下重建分析输入。
#
# 【Purpose / 目的】
# 从GEO下载supplementary files，或确认本地输入目录。
# Download GEO supplementary files or validate a local input directory.
#
# 【为什么优先使用processed raw counts？ / Why processed raw counts first?】
# 对学习和常规下游scRNA-seq分析，作者提供的filtered/raw count matrix
# 通常足够；不必一开始就从FASTQ重新跑Cell Ranger。
# For downstream learning/analysis, author-provided count matrices are
# usually sufficient; re-running Cell Ranger from FASTQ is not mandatory.
#
# 【注意 / Important】
# “processed raw counts”仍然可以是整数UMI counts。
# “Processed” in GEO does not necessarily mean normalized expression.
#
# 【Output / 输出】
# 输入文件清单和M02_input_location.rds。
# An input-file manifest and M02_input_location.rds.
#
# 【Check after running / 运行后检查】
# 检查文件数、文件扩展名和目录是否合理。
# Check file count, extensions and paths before reading matrices.
# ============================================================

if (RUN_M02_ACQUIRE_DATA) {

  MODULE <- "M02_acquire_data"
  MDIR <- module_dir(MODULE)

  if (DATA_SOURCE == "GEO") {

    if (!nzchar(GSE_ID)) {
      safe_stop_module(
        MODULE,
        "DATA_SOURCE='GEO'时必须设置GSE_ID。"
      )
    }

    if (!requireNamespace("GEOquery", quietly = TRUE)) {
      safe_stop_module(
        MODULE,
        "缺少GEOquery。请先运行M00安装。"
      )
    }

    gse_dir <- file.path(
      DIR_RAW,
      GSE_ID
    )

    dir.create(
      gse_dir,
      recursive = TRUE,
      showWarnings = FALSE
    )

    if (DOWNLOAD_GEO_DATA) {

      message(
        "Downloading GEO supplementary files: ",
        GSE_ID
      )

      GEOquery::getGEOSuppFiles(
        GEO = GSE_ID,
        baseDir = DIR_RAW,
        makeDirectory = TRUE,
        fetch_files = TRUE
      )
    }


    downloaded <- list.files(
      gse_dir,
      recursive = TRUE,
      full.names = TRUE
    )

    archives <- downloaded[
      grepl(
        "\\.(tar|tar\\.gz|tgz)$",
        downloaded,
        ignore.case = TRUE
      )
    ]

    if (length(archives) > 0L) {

      extracted_dir <- file.path(
        gse_dir,
        "extracted"
      )

      dir.create(
        extracted_dir,
        recursive = TRUE,
        showWarnings = FALSE
      )

      existing_extract <- list.files(
        extracted_dir,
        recursive = TRUE,
        full.names = TRUE
      )

      if (
        FORCE_REEXTRACT ||
        length(existing_extract) == 0L
      ) {

        for (archive in archives) {

          message(
            "Extracting: ",
            archive
          )

          utils::untar(
            archive,
            exdir = extracted_dir
          )
        }
      }

      LOCAL_INPUT_DIR <- extracted_dir

    } else {

      # 有些GEO数据不是tar，而是直接多个补充文件。
      LOCAL_INPUT_DIR <- gse_dir
    }

  } else if (DATA_SOURCE == "LOCAL") {

    if (!dir.exists(LOCAL_INPUT_DIR)) {
      safe_stop_module(
        MODULE,
        paste0(
          "本地输入目录不存在 / Local input directory not found:\n",
          LOCAL_INPUT_DIR,
          "\n\n请把输入文件放到该目录，或修改第0节INPUT_DIR。\n",
          "Place the input files there, or edit INPUT_DIR in Section 0."
        )
      )
    }

  } else {
    safe_stop_module(
      MODULE,
      "DATA_SOURCE只能是'GEO'或'LOCAL'。"
    )
  }


  all_files <- list.files(
    LOCAL_INPUT_DIR,
    recursive = TRUE,
    full.names = TRUE
  )

  writeLines(
    all_files,
    file.path(
      MDIR,
      "tables",
      "input_file_manifest.txt"
    )
  )

  cat(
    "\nInput directory:\n",
    LOCAL_INPUT_DIR,
    "\nFiles found: ",
    length(all_files),
    "\n"
  )

  save_checkpoint(
    list(
      input_dir = LOCAL_INPUT_DIR,
      files = all_files
    ),
    "M02_input_location"
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M03. INPUT AUDIT / 输入审计
#
# 【Module scope / 模块范围】
# The audit tests whether file organization is internally consistent before large
# matrices are loaded. Passing this module confirms structural readability, not
# biological correctness or sample-label validity.
# 本模块在加载大型矩阵前检查文件组织是否内部一致。通过审计仅表示结构上可以读取，
# 不代表生物学内容正确，也不代表样本标签已经得到验证。
#
# 【Stop condition / 停止条件】
# Missing matrix triplets, ambiguous sample-to-file mappings or unexpected feature
# types must be resolved from source documentation before proceeding. Automatic
# guessing would make downstream provenance unreliable.
# 若矩阵三件套缺失、样本与文件对应不唯一或feature type异常，应依据源数据说明解决
# 后再继续；自动猜测会破坏下游分析的可追溯性。
#
# 【Purpose / 目的】
# 在真正创建Seurat对象之前先检查文件结构。
# Inspect the file structure before constructing Seurat objects.
#
# 【主要检查 / Main checks】
# 1. matrix.mtx、features.tsv、barcodes.tsv是否成套存在。
#    Whether MTX/features/barcodes occur as valid triplets.
# 2. features.tsv是否包含第三列feature type。
#    Whether a feature-type column is available.
# 3. 是否存在Gene Expression / Antibody Capture / HTO等feature。
#    Whether RNA and/or antibody/HTO features are present.
#
# 【为什么单独设这一模块 / Why a separate audit?】
# 很多GEO数据读取失败不是分析算法问题，而是文件命名/结构不同。
# Many GEO failures arise from file organization rather than biology.
# 先审计可以在最早阶段发现问题。
# Auditing catches those problems before QC or clustering.
#
# 【Output / 输出】
#   input_audit.csv
#   preflight_validation.csv
#   input_preflight_validation.png
#   input_feature_type_composition.png（仅在feature-type计数可用时）
#   M03_input_audit.rds
#
# 【Check after running / 运行后检查】
# 每个sample是否只有一套正确matrix？HTO是否真的存在？
# Does each sample map to the expected matrix, and is HTO really present?
# ============================================================

if (RUN_M03_INPUT_AUDIT) {

  MODULE <- "M03_input_audit"
  MDIR <- module_dir(MODULE)

  if (!exists("sample_metadata")) {
    sample_metadata <- load_checkpoint(
      "M01_sample_metadata"
    )
  }

  loc <- load_checkpoint(
    "M02_input_location"
  )

  LOCAL_INPUT_DIR <- loc$input_dir


  audit_rows <- list()


  if (INPUT_FORMAT == "10X_MTX") {

    matrix_files <- list.files(
      LOCAL_INPUT_DIR,
      pattern = "matrix\\.mtx(\\.gz)?$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (length(matrix_files) == 0L) {
      safe_stop_module(
        MODULE,
        paste0(
          "没有找到matrix.mtx(.gz)：\n",
          LOCAL_INPUT_DIR
        )
      )
    }


    for (matrix_file in matrix_files) {

      stem <- sub(
        "matrix\\.mtx(\\.gz)?$",
        "",
        matrix_file,
        ignore.case = TRUE
      )

      feature_file <- first_existing(
        c(
          paste0(
            stem,
            "features.tsv.gz"
          ),
          paste0(
            stem,
            "features.tsv"
          ),
          paste0(
            stem,
            "genes.tsv.gz"
          ),
          paste0(
            stem,
            "genes.tsv"
          )
        )
      )

      barcode_file <- first_existing(
        c(
          paste0(
            stem,
            "barcodes.tsv.gz"
          ),
          paste0(
            stem,
            "barcodes.tsv"
          )
        )
      )

      gsm <- extract_gsm(
        basename(matrix_file)
      )

      fallback_sample <- clean_file_stem(
        matrix_file
      )

      if (!nzchar(fallback_sample)) {

        parent_path <- dirname(matrix_file)
        generic_10x_dirs <- c(
          "filtered_feature_bc_matrix",
          "raw_feature_bc_matrix",
          "outs"
        )

        repeat {

          parent_name <- basename(parent_path)

          if (
            nzchar(parent_name) &&
            !parent_name %in% generic_10x_dirs
          ) {
            break
          }

          next_parent <- dirname(parent_path)

          if (identical(next_parent, parent_path)) {
            break
          }

          parent_path <- next_parent
        }

        fallback_sample <- parent_name
      }

      sample_id <- if (
        !is.na(gsm)
      ) {
        gsm
      } else {
        fallback_sample
      }


      n_features <- NA_integer_
      feature_types <- ""
      n_gene_expression <- NA_integer_
      n_antibody <- NA_integer_
      hto_like_features <- ""


      if (
        !is.na(feature_file) &&
        file.exists(feature_file)
      ) {

        ft <- data.table::fread(
          feature_file,
          header = FALSE,
          data.table = FALSE
        )

        n_features <- nrow(ft)

        if (ncol(ft) >= 3L) {

          feature_types <- paste(
            unique(
              as.character(
                ft[[3L]]
              )
            ),
            collapse = "; "
          )

          n_gene_expression <- sum(
            grepl(
              "Gene Expression",
              ft[[3L]],
              ignore.case = TRUE
            )
          )

          # Known 10x feature-barcode labels include "Antibody Capture".
          # tissue hashtag features as "Custom".
          #
          # 已知10x feature barcode常标为Antibody Capture；
          is_gene_expression <- grepl(
            "^Gene Expression$",
            ft[[3L]],
            ignore.case = TRUE
          )

          feature_type_values <- as.character(ft[[3L]])

          is_known_hto <- feature_type_values %in% HTO_FEATURE_TYPES |
            grepl(
              "HTO|Hashtag",
              feature_type_values,
              ignore.case = TRUE
            )

          # v6 safety rule:
          # non-GEX features are NOT automatically HTO unless explicitly enabled.
          is_hto_candidate <- is_known_hto |
            (HAS_HTO & ALLOW_NON_GEX_AS_HTO & !is_gene_expression)

          n_antibody <- sum(
            is_hto_candidate
          )

          if (any(is_hto_candidate)) {

            hto_names <- unique(
              as.character(
                ft[
                  is_hto_candidate,
                  min(
                    2L,
                    ncol(ft)
                  )
                ]
              )
            )

            hto_like_features <- paste(
              hto_names,
              collapse = "; "
            )
          }

        } else {
          feature_types <- "feature_type_column_absent"
        }
      }


      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        sample_id = sample_id,
        matrix_file = matrix_file,
        feature_file = feature_file,
        barcode_file = barcode_file,
        n_features = n_features,
        feature_types = feature_types,
        n_gene_expression =
          n_gene_expression,
        n_antibody_or_hto =
          n_antibody,
        hto_like_features =
          hto_like_features,
        stringsAsFactors = FALSE
      )
    }


  } else if (INPUT_FORMAT == "10X_H5") {

    h5_files <- list.files(
      LOCAL_INPUT_DIR,
      pattern = "\\.h5$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )

    if (
      nzchar(INPUT_H5_FILE) &&
      file.exists(INPUT_H5_FILE)
    ) {
      h5_files <- unique(
        c(
          INPUT_H5_FILE,
          h5_files
        )
      )
    }

    if (length(h5_files) == 0L) {
      safe_stop_module(
        MODULE,
        "没有找到10x H5文件。"
      )
    }

    audit_rows <- lapply(
      h5_files,
      function(x) {

        data.frame(
          sample_id = ifelse(
            !is.na(
              extract_gsm(
                basename(x)
              )
            ),
            extract_gsm(
              basename(x)
            ),
            tools::file_path_sans_ext(
              basename(x)
            )
          ),
          h5_file = x,
          stringsAsFactors = FALSE
        )
      }
    )


  } else if (
    INPUT_FORMAT %in%
      c(
        "SEURAT_RDS",
        "SCE_RDS"
      )
  ) {

    if (
      !nzchar(INPUT_RDS_FILE) ||
      !file.exists(INPUT_RDS_FILE)
    ) {
      safe_stop_module(
        MODULE,
        "INPUT_FORMAT为RDS时，请设置存在的INPUT_RDS_FILE。"
      )
    }

    audit_rows <- list(
      data.frame(
        sample_id = "RDS_object",
        rds_file = INPUT_RDS_FILE,
        stringsAsFactors = FALSE
      )
    )


  } else {
    safe_stop_module(
      MODULE,
      "不支持的INPUT_FORMAT。"
    )
  }


  input_audit <- data.table::rbindlist(
    audit_rows,
    fill = TRUE
  )

  input_audit <- as.data.frame(
    input_audit,
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    input_audit,
    file.path(
      MDIR,
      "tables",
      "input_audit.csv"
    )
  )

  print(input_audit)

  # 对HTO模式进行提示，但不在这里武断停止：
  # 某些HTO counts可能位于独立文件而不是RNA features.tsv。
  if (
    MULTIPLEXING_MODE == "HTO" &&
    "n_antibody_or_hto" %in%
      colnames(input_audit)
  ) {

    if (
      all(
        is.na(
          input_audit$n_antibody_or_hto
        ) |
        input_audit$n_antibody_or_hto == 0L
      )
    ) {
      warning(
        "当前MTX features.tsv没有检测到Antibody/HTO feature。\n",
        "这不一定证明HTO不存在；可能存放在独立feature-barcode文件中。\n",
        "在进入M05前必须确认组织标签来源。"
      )
    }
  }

  # ---------- v6 strict preflight validation ----------
  preflight <- data.frame(
    check = character(0),
    status = character(0),
    message = character(0),
    stringsAsFactors = FALSE
  )

  add_check <- function(check, ok, message_ok, message_fail) {
    preflight <<- rbind(
      preflight,
      data.frame(
        check = check,
        status = if (isTRUE(ok)) "PASS" else "FAIL",
        message = if (isTRUE(ok)) message_ok else message_fail,
        stringsAsFactors = FALSE
      )
    )
  }

  if (INPUT_FORMAT == "10X_MTX") {
    add_check(
      "unique_sample_id",
      !anyDuplicated(input_audit$sample_id),
      "Each matrix maps to a unique sample_id.",
      "Duplicated sample_id detected in input audit."
    )

    add_check(
      "triplet_files_exist",
      all(
        !is.na(input_audit$matrix_file) & file.exists(input_audit$matrix_file) &
        !is.na(input_audit$feature_file) & file.exists(input_audit$feature_file) &
        !is.na(input_audit$barcode_file) & file.exists(input_audit$barcode_file)
      ),
      "All matrix/features/barcodes files exist.",
      "At least one MTX/features/barcodes triplet is incomplete."
    )

    input_ids <- sort(unique(as.character(input_audit$sample_id)))
    metadata_ids <- sort(unique(as.character(sample_metadata[[SAMPLE_ID_COL]])))

    add_check(
      "input_metadata_match",
      identical(input_ids, metadata_ids),
      "Input sample IDs exactly match sample metadata.",
      paste0(
        "Input/metadata sample mismatch. Missing in metadata: ",
        paste(setdiff(input_ids, metadata_ids), collapse = ", "),
        "; missing in input: ",
        paste(setdiff(metadata_ids, input_ids), collapse = ", ")
      )
    )

    add_check(
      "gene_expression_present",
      all(is.na(input_audit$n_gene_expression) | input_audit$n_gene_expression > 0L),
      "Gene Expression features are present.",
      "At least one sample has no Gene Expression features."
    )

    if (MULTIPLEXING_MODE == "HTO") {
      add_check(
        "hto_features_present",
        all(!is.na(input_audit$n_antibody_or_hto) & input_audit$n_antibody_or_hto > 0L),
        "HTO candidate features are present in every sample.",
        paste0(
          "HTO mode is enabled but HTO features were not detected in every sample. ",
          "Check HTO_FEATURE_TYPES / profile settings before M04-M05."
        )
      )
    }
  }

  add_check(
    "biological_replicate_complete",
    all(
      !is.na(sample_metadata[[BIOLOGICAL_REPLICATE_COL]]) &
      trimws(as.character(sample_metadata[[BIOLOGICAL_REPLICATE_COL]])) != ""
    ),
    "Biological replicate IDs are complete.",
    "Missing biological replicate IDs detected."
  )

  add_check(
    "condition_complete",
    all(
      !is.na(sample_metadata[[CONDITION_COL]]) &
      trimws(as.character(sample_metadata[[CONDITION_COL]])) != ""
    ),
    "Condition labels are complete.",
    "Missing condition labels detected."
  )

  add_check(
    "control_case_present",
    all(c(CONTROL_GROUP, CASE_GROUP) %in% unique(sample_metadata[[CONDITION_COL]])),
    "Control and case groups are both present.",
    "CONTROL_GROUP and/or CASE_GROUP are absent from metadata."
  )

  data.table::fwrite(
    preflight,
    file.path(MDIR, "tables", "preflight_validation.csv")
  )


  # ==========================================================
  # INPUT PREFLIGHT FIGURE / 输入预检图
  #
  # 该图直接使用preflight_validation.csv中的检查结果，不从文件名或
  # 示意图推断数据内容。可见文字统一使用英文，避免不同PNG设备缺少
  # 中文字体时出现方框或乱码；源代码注释仍保持中英双语。
  #
  # This figure is generated directly from the recorded preflight checks.
  # Visible plot text is English-only so that PNG rendering does not depend
  # on the availability of a Chinese font on the current system.
  # ==========================================================

  preflight_plot_data <- preflight |>
    dplyr::mutate(
      check_label = gsub(
        "_",
        " ",
        .data$check,
        fixed = TRUE
      ),
      check_label = factor(
        .data$check_label,
        levels = rev(
          unique(
            .data$check_label
          )
        )
      ),
      status = factor(
        .data$status,
        levels = c(
          "PASS",
          "FAIL"
        )
      )
    )

  p_input_preflight <- ggplot2::ggplot(
    preflight_plot_data,
    ggplot2::aes(
      x = 1,
      y = .data$check_label,
      fill = .data$status
    )
  ) +
    ggplot2::geom_tile(
      width = 0.82,
      height = 0.72,
      colour = "white",
      linewidth = 0.7
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = .data$status
      ),
      colour = "white",
      fontface = "bold",
      size = 3.8
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "PASS" = "#2E8B57",
        "FAIL" = "#C23B22"
      ),
      drop = FALSE,
      name = NULL
    ) +
    ggplot2::labs(
      title = "Input preflight validation",
      subtitle = paste0(
        "Input format: ",
        INPUT_FORMAT,
        "; checks are evaluated before expression matrices enter analysis"
      ),
      x = NULL,
      y = NULL,
      caption = "A PASS result confirms structural consistency, not biological correctness."
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      legend.position = "none",
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )

  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "input_preflight_validation.png"
    ),
    p_input_preflight,
    width = 8.5,
    height = max(
      4.8,
      0.48 * nrow(preflight_plot_data) + 2.2
    ),
    dpi = 300
  )


  # 仅当M03已从10x features.tsv获得feature-type计数时绘制组成图。
  # For other input formats, M04 will report the actual object dimensions;
  # M03 does not load a large RDS merely to imitate a feature-composition plot.
  if (
    INPUT_FORMAT == "10X_MTX" &&
    all(
      c(
        "n_features",
        "n_gene_expression",
        "n_antibody_or_hto"
      ) %in% colnames(input_audit)
    )
  ) {

    feature_composition <- input_audit |>
      dplyr::transmute(
        sample_id = as.character(
          .data$sample_id
        ),
        `Gene Expression` = as.numeric(
          .data$n_gene_expression
        ),
        `Antibody or HTO` = as.numeric(
          .data$n_antibody_or_hto
        ),
        Other = pmax(
          as.numeric(
            .data$n_features
          ) -
            dplyr::coalesce(
              as.numeric(
                .data$n_gene_expression
              ),
              0
            ) -
            dplyr::coalesce(
              as.numeric(
                .data$n_antibody_or_hto
              ),
              0
            ),
          0
        )
      ) |>
      tidyr::pivot_longer(
        cols = -sample_id,
        names_to = "feature_type",
        values_to = "n_features"
      ) |>
      dplyr::filter(
        is.finite(
          .data$n_features
        ),
        .data$n_features >= 0
      )

    if (
      nrow(feature_composition) > 0L &&
      sum(
        feature_composition$n_features,
        na.rm = TRUE
      ) > 0
    ) {

      p_feature_composition <- ggplot2::ggplot(
        feature_composition,
        ggplot2::aes(
          x = .data$sample_id,
          y = .data$n_features,
          fill = .data$feature_type
        )
      ) +
        ggplot2::geom_col(
          position = "fill",
          width = 0.72,
          colour = "white",
          linewidth = 0.25
        ) +
        ggplot2::scale_fill_manual(
          values = c(
            "Gene Expression" = "#3C78D8",
            "Antibody or HTO" = "#F6A13A",
            "Other" = "#4FAF72"
          ),
          name = "Feature type"
        ) +
        ggplot2::scale_y_continuous(
          labels = function(x) {
            paste0(
              round(
                100 * x
              ),
              "%"
            )
          },
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.04
            )
          )
        ) +
        ggplot2::labs(
          title = "Input feature-type composition",
          subtitle = "Composition is calculated from the third column of the 10x feature table",
          x = NULL,
          y = "Feature proportion",
          caption = "Non-gene-expression features are counted as HTO only under the configured HTO rules."
        ) +
        ggplot2::theme_minimal(
          base_size = 11
        ) +
        ggplot2::theme(
          panel.grid.major.x = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          ),
          legend.position = "bottom"
        )

      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "input_feature_type_composition.png"
        ),
        p_feature_composition,
        width = max(
          7.5,
          0.62 *
            dplyr::n_distinct(
              feature_composition$sample_id
            ) +
            3.5
        ),
        height = 6,
        dpi = 300
      )
    }
  }

  if (any(preflight$status == "FAIL")) {
    failed <- preflight$message[preflight$status == "FAIL"]
    safe_stop_module(
      MODULE,
      paste0(
        "M03 preflight validation failed:\n- ",
        paste(failed, collapse = "\n- ")
      )
    )
  }

  save_checkpoint(
    input_audit,
    "M03_input_audit"
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M04. READ EXPRESSION INPUT / 读取表达矩阵
#
# 【Module scope / 模块范围】
# Supported inputs are converted to a consistent Seurat representation while raw
# count information and sample provenance are preserved. Conversion harmonizes
# object structure; it does not normalize expression or remove cells.
# 本模块在保留原始计数和样本来源信息的前提下，将不同输入统一为Seurat对象结构。
# 结构统一不等同于表达标准化，也不会在此阶段删除细胞。
#
# 【Quality gate / 质量门槛】
# Gene-by-cell dimensions, assay names, barcode uniqueness and raw-count
# availability must be plausible for every sample. Unexpected fractional values
# in the designated raw-count layer require investigation before count-based DE.
# 每个样本的基因×细胞维度、assay名称、barcode唯一性及原始计数可用性都应合理。
# 若指定原始计数layer中出现异常小数值，必须在基于计数的差异分析前查明原因。
#
# 【Purpose / 目的】
# 把不同输入格式统一转换为Seurat object/list。
# Convert supported input formats into a consistent Seurat object/list.
#
# 【核心对象 / Core object】
# RNA assay:
#   rows = genes
#   columns = cells
#   values = raw UMI/read counts
#
# 可选HTO assay / optional HTO assay:
#   rows = hashtag/antibody features
#   columns = the same cell barcodes
#
# 【为什么保留raw counts / Why retain raw counts?】
# 后面的normalization会产生用于PCA/UMAP的表达值，
# 但pseudobulk edgeR必须回到原始整数counts。
# Normalization creates values for exploratory analyses, while
# pseudobulk edgeR must use raw integer counts.
#
# 【Output / 输出】
#   object_summary.csv
#   sample_summary.csv（当sample_id metadata可用时）
#   input_object_dimensions.png
#   input_cells_by_sample.png（当sample_id metadata可用时）
#   M04_seurat_list.rds
#
# 【Check after running / 运行后检查】
# 每个sample的细胞数、基因数、assay是否符合预期。
# Check cell number, feature number and assay names for every sample.
# ============================================================

if (RUN_M04_READ_INPUT) {

  MODULE <- "M04_read_input"
  MDIR <- module_dir(MODULE)

  if (!exists("sample_metadata")) {
    sample_metadata <- load_checkpoint(
      "M01_sample_metadata"
    )
  }

  input_audit <- load_checkpoint(
    "M03_input_audit"
  )


  add_sample_level_metadata <- function(
      object,
      sample_id
  ) {

    idx <- match(
      sample_id,
      sample_metadata[[SAMPLE_ID_COL]]
    )

    if (is.na(idx)) {
      stop(
        "sample_id未在sample_metadata中找到：",
        sample_id
      )
    }

    object$sample_id <- sample_id

    for (
      col_name in
      colnames(sample_metadata)
    ) {

      object[[col_name]] <- sample_metadata[[col_name]][idx]
    }

    object
  }


  read_mtx_with_feature_types <- function(
      matrix_file,
      feature_file,
      barcode_file,
      sample_id
  ) {

    raw_matrix <- Seurat::ReadMtx(
      mtx = matrix_file,
      cells = barcode_file,
      features = feature_file,
      cell.column = 1,
      feature.column = 2,
      unique.features = TRUE,
      strip.suffix = FALSE
    )

    feature_table <- data.table::fread(
      feature_file,
      header = FALSE,
      data.table = FALSE
    )

    if (nrow(feature_table) != nrow(raw_matrix)) {
      stop(
        "feature行数与matrix行数不一致：",
        sample_id
      )
    }


    if (ncol(feature_table) >= 3L) {

      feature_type <- as.character(
        feature_table[[3L]]
      )

      gene_rows <- grepl(
        "^Gene Expression$",
        feature_type,
        ignore.case = TRUE
      )

      known_hto_rows <- feature_type %in% HTO_FEATURE_TYPES |
        grepl(
          "HTO|Hashtag",
          feature_type,
          ignore.case = TRUE
        )

      if (!any(gene_rows)) {
        # If no standard "Gene Expression" label exists, exclude known
        # feature-barcode rows and conservatively treat the rest as RNA.
        # 若不存在标准Gene Expression标签，则排除已知HTO后把其余feature作为RNA候选。
        gene_rows <- !known_hto_rows
      }

      # Some GEO submissions store HTOs as feature_type="Custom".
      # When HAS_HTO=TRUE, non-Gene-Expression rows are therefore retained
      # as an HTO assay rather than silently discarded.
      #
      # 部分GEO数据将HTO标为Custom。
      # 当HAS_HTO=TRUE时，非Gene Expression行会保留到HTO assay中。
      antibody_rows <- (
        known_hto_rows |
          (HAS_HTO & ALLOW_NON_GEX_AS_HTO & !gene_rows)
      )

    } else {

      gene_rows <- rep(
        TRUE,
        nrow(raw_matrix)
      )

      antibody_rows <- rep(
        FALSE,
        nrow(raw_matrix)
      )
    }


    rna_counts <- raw_matrix[
      gene_rows,
      ,
      drop = FALSE
    ]

    object <- Seurat::CreateSeuratObject(
      counts = rna_counts,
      project = sample_id,
      min.cells = 0,
      min.features = 0
    )

    object <- Seurat::RenameCells(
      object,
      add.cell.id = sample_id
    )


    if (any(antibody_rows)) {

      hto_counts <- raw_matrix[
        antibody_rows,
        ,
        drop = FALSE
      ]

      # RenameCells已经给RNA细胞加前缀；
      # HTO矩阵列名必须同步。
      colnames(hto_counts) <- paste0(
        sample_id,
        "_",
        colnames(hto_counts)
      )

      object[["HTO"]] <- SeuratObject::CreateAssayObject(
        counts = hto_counts
      )
    }


    object <- add_sample_level_metadata(
      object,
      sample_id
    )

    object
  }


  seurat_list <- list()


  if (INPUT_FORMAT == "10X_MTX") {

    for (
      i in seq_len(
        nrow(input_audit)
      )
    ) {

      row <- input_audit[
        i,
        ,
        drop = FALSE
      ]

      sample_id <- as.character(
        row$sample_id
      )

      if (
        !sample_id %in%
          sample_metadata[[SAMPLE_ID_COL]]
      ) {

        # 对自己的实验，文件名可能无法自动推断sample_id。
        # 优先尝试按sample_metadata中的input_path匹配。
        warning(
          "自动推断的sample_id不在sample_metadata：",
          sample_id,
          "\n请检查文件名或在sample_metadata中使用一致sample_id。"
        )
      }

      obj <- read_mtx_with_feature_types(
        matrix_file =
          row$matrix_file,
        feature_file =
          row$feature_file,
        barcode_file =
          row$barcode_file,
        sample_id =
          sample_id
      )

      seurat_list[[sample_id]] <- obj
    }


  } else if (INPUT_FORMAT == "10X_H5") {

    for (
      i in seq_len(
        nrow(input_audit)
      )
    ) {

      h5_file <- input_audit$h5_file[i]
      sample_id <- input_audit$sample_id[i]

      h5 <- Seurat::Read10X_h5(
        h5_file,
        use.names = TRUE,
        unique.features = TRUE
      )

      if (is.list(h5)) {

        rna_name <- names(h5)[
          grepl(
            "Gene Expression",
            names(h5),
            ignore.case = TRUE
          )
        ][1]

        if (is.na(rna_name)) {
          rna_name <- names(h5)[1]
        }

        obj <- CreateSeuratObject(
          counts = h5[[rna_name]],
          project = sample_id,
          min.cells = 0,
          min.features = 0
        )

        antibody_name <- names(h5)[
          grepl(
            "Antibody|HTO",
            names(h5),
            ignore.case = TRUE
          )
        ][1]

        if (!is.na(antibody_name)) {
          obj[["HTO"]] <-
            CreateAssayObject(
              counts = h5[[antibody_name]])
        }

      } else {

        obj <- CreateSeuratObject(
          counts = h5,
          project = sample_id,
          min.cells = 0,
          min.features = 0
        )
      }

      obj <- RenameCells(
        obj,
        add.cell.id = sample_id
      )

      obj <- add_sample_level_metadata(
        obj,
        sample_id
      )

      seurat_list[[sample_id]] <- obj
    }


  } else if (
    INPUT_FORMAT == "SEURAT_RDS"
  ) {

    obj <- readRDS(
      INPUT_RDS_FILE
    )

    if (!inherits(obj, "Seurat")) {
      safe_stop_module(
        MODULE,
        "INPUT_RDS_FILE不是Seurat对象。"
      )
    }

    seurat_list <- list(
      RDS_object = obj
    )


  } else if (
    INPUT_FORMAT == "SCE_RDS"
  ) {

    sce <- readRDS(
      INPUT_RDS_FILE
    )

    if (
      !inherits(
        sce,
        "SingleCellExperiment"
      )
    ) {
      safe_stop_module(
        MODULE,
        "INPUT_RDS_FILE不是SingleCellExperiment对象。"
      )
    }

    obj <- as.Seurat(
      sce,
      counts = "counts"
    )

    seurat_list <- list(
      RDS_object = obj
    )
  }


  save_checkpoint(
    seurat_list,
    "M04_seurat_list"
  )


  # ============================================================
  # OBJECT-LEVEL SUMMARY / 输入对象层面汇总
  #
  # ncol()/nrow()在部分Seurat/SeuratObject版本中可能返回double，
  # 因此FUN.VALUE使用numeric(1)，避免integer/double类型不一致。
  #
  # For prepared SEURAT_RDS input, one RDS object may contain
  # multiple biological samples. Therefore object_id describes the
  # input object and should not automatically be interpreted as sample_id.
  # ============================================================

  object_summary <- data.frame(
    object_id = names(seurat_list),

    n_cells = vapply(
      seurat_list,
      function(x) as.numeric(ncol(x)),
      FUN.VALUE = numeric(1)
    ),

    n_genes = vapply(
      seurat_list,
      function(x) as.numeric(nrow(x)),
      FUN.VALUE = numeric(1)
    ),

    assays = vapply(
      seurat_list,
      function(x) {
        paste(
          SeuratObject::Assays(x),
          collapse = "; "
        )
      },
      FUN.VALUE = character(1)
    ),

    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    object_summary,
    file.path(
      MDIR,
      "tables",
      "object_summary.csv"
    )
  )

  print(
    object_summary
  )


  # ============================================================
  # TRUE SAMPLE-LEVEL SUMMARY / 真实样本层面汇总
  #
  # prepared Seurat RDS内部可能已经包含多个真实生物学样本。
  # 因此进一步读取cell-level metadata中的SAMPLE_ID_COL，
  # 汇总真正的sample组成和每个sample的细胞数。
  # ============================================================

  sample_summary <- NULL

  sample_summary_list <- lapply(
    seurat_list,
    function(x) {

      metadata <- x[[]]

      if (!(SAMPLE_ID_COL %in% colnames(metadata))) {
        return(NULL)
      }

      sample_values <- as.character(
        metadata[, SAMPLE_ID_COL]
      )

      sample_counts <- table(
        sample_values,
        useNA = "ifany"
      )

      data.frame(
        sample_id = names(sample_counts),
        n_cells = as.numeric(sample_counts),
        stringsAsFactors = FALSE
      )
    }
  )

  sample_summary_list <- Filter(
    Negate(is.null),
    sample_summary_list
  )

  if (length(sample_summary_list) > 0L) {

    sample_summary <- do.call(
      rbind,
      sample_summary_list
    )

    rownames(sample_summary) <- NULL

    data.table::fwrite(
      sample_summary,
      file.path(
        MDIR,
        "tables",
        "sample_summary.csv"
      )
    )

    cat(
      "\nTrue sample-level summary / 真实sample层面汇总:\n"
    )

    print(
      sample_summary
    )

  } else {

    warning(
      "No sample_id metadata found / 没有找到sample_id metadata。"
    )
  }


  # ==========================================================
  # INPUT OBJECT FIGURES / 输入对象图片
  #
  # M04尚未进行normalization、PCA或UMAP，因此这里不绘制“看起来像UMAP”
  # 的示意散点图。以下图片只展示刚读入对象的真实维度和真实样本细胞数，
  # 从而保持模块边界与结果可追溯性。
  #
  # M04 has not normalized the data or run PCA/UMAP. Therefore these plots
  # show the measured object dimensions and sample cell counts instead of a
  # synthetic embedding that could be mistaken for an analysis result.
  # ==========================================================

  object_dimension_plot_data <- object_summary |>
    dplyr::select(
      object_id,
      n_cells,
      n_genes
    ) |>
    tidyr::pivot_longer(
      cols = c(
        n_cells,
        n_genes
      ),
      names_to = "metric",
      values_to = "value"
    ) |>
    dplyr::mutate(
      metric = factor(
        .data$metric,
        levels = c(
          "n_cells",
          "n_genes"
        ),
        labels = c(
          "Cells",
          "Genes"
        )
      )
    )

  p_object_dimensions <- ggplot2::ggplot(
    object_dimension_plot_data,
    ggplot2::aes(
      x = stats::reorder(
        .data$object_id,
        .data$value
      ),
      y = .data$value,
      fill = .data$metric
    )
  ) +
    ggplot2::geom_col(
      width = 0.68,
      show.legend = FALSE
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        label = format(
          .data$value,
          big.mark = ",",
          scientific = FALSE,
          trim = TRUE
        )
      ),
      hjust = -0.12,
      size = 3.2
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(
        metric
      ),
      scales = "free_y",
      ncol = 1
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Cells" = "#3C78D8",
        "Genes" = "#54A76B"
      )
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.18
        )
      )
    ) +
    ggplot2::labs(
      title = "Expression-input object dimensions",
      subtitle = paste0(
        "Objects were read as ",
        INPUT_FORMAT,
        "; no normalization or cell filtering has been applied"
      ),
      x = NULL,
      y = "Count",
      caption = "Object identifiers describe input containers and are not automatically biological sample IDs."
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 30,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      )
    )

  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "input_object_dimensions.png"
    ),
    p_object_dimensions,
    width = max(
      7.5,
      0.62 * nrow(object_summary) + 4
    ),
    height = 7.5,
    dpi = 300
  )


  if (
    !is.null(sample_summary) &&
    nrow(sample_summary) > 0L
  ) {

    sample_plot_data <- sample_summary |>
      dplyr::group_by(
        sample_id
      ) |>
      dplyr::summarise(
        n_cells = sum(
          .data$n_cells
        ),
        .groups = "drop"
      ) |>
      dplyr::left_join(
        sample_metadata |>
          dplyr::transmute(
            sample_id = as.character(
              .data[[SAMPLE_ID_COL]]
            ),
            condition = as.character(
              .data[[CONDITION_COL]]
            )
          ) |>
          dplyr::distinct(),
        by = "sample_id"
      ) |>
      dplyr::mutate(
        condition = dplyr::coalesce(
          .data$condition,
          "Unspecified"
        )
      )

    p_cells_by_sample <- ggplot2::ggplot(
      sample_plot_data,
      ggplot2::aes(
        x = stats::reorder(
          .data$sample_id,
          .data$n_cells
        ),
        y = .data$n_cells,
        fill = .data$condition
      )
    ) +
      ggplot2::geom_col(
        width = 0.70
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = format(
            .data$n_cells,
            big.mark = ",",
            scientific = FALSE,
            trim = TRUE
          )
        ),
        hjust = -0.12,
        size = 3.1
      ) +
      ggplot2::coord_flip(
        clip = "off"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(
          mult = c(
            0,
            0.18
          )
        )
      ) +
      ggplot2::labs(
        title = "Cells recovered from each biological sample",
        subtitle = "Counts are read from cell-level sample metadata before QC and doublet removal",
        x = NULL,
        y = "Number of cells",
        fill = "Condition",
        caption = "Large differences should be checked against library loading and sequencing depth."
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "input_cells_by_sample.png"
      ),
      p_cells_by_sample,
      width = 8.5,
      height = max(
        5.5,
        0.38 * nrow(sample_plot_data) + 3.2
      ),
      dpi = 300
    )
  }


  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M05. OPTIONAL HTO/CELL HASHING DEMULTIPLEXING / 可选HTO/Cell Hashing解复用
#
# 【Module scope / 模块范围】
# This optional module resolves pooled identities from HTO counts or imports an
# externally validated cell map. It is not needed when each library already
# represents one known specimen and tissue.
# 本可选模块利用HTO计数拆分混合来源，或导入已经验证的细胞级映射。若每个文库本身
# 已对应一个明确的标本和组织，则不需要运行该步骤。
#
# 【Interpretation boundary / 解读边界】
# HTO singlet/doublet/negative labels describe hashtag assignment, not RNA-based
# transcriptomic doublets. HTO identities become biological sample or tissue labels
# only after HTO_LABEL_MAP and HTO_LABEL_ROLE are verified from the experiment.
# HTO的singlet/doublet/negative描述hashtag归属，并不等同于基于RNA表达识别的
# transcriptomic doublet。只有结合实验记录验证HTO_LABEL_MAP和HTO_LABEL_ROLE后，
# HTO identity才能解释为生物学样本或组织标签。
#
# 【何时需要 / When is this needed?】
# 只有实验在建库前把多个样本/组织pool在一起，并使用HTO/Cell Hashing
# 标记来源时才需要。
# Use this only when cells from multiple samples/tissues were pooled
# and tagged by hashtags before sequencing.
#
# 普通“每只小鼠一个独立Liver library”的实验应设置：
# For ordinary one-liver-library-per-animal experiments:
#   MULTIPLEXING_MODE <- "NONE"
#
# 【步骤 / Steps】
# HTO counts → CLR normalization → HTODemux →
# Negative / Singlet / Doublet classification.
#
# 【注意 / Important】
# HTO doublet和scDblFinder doublet不是完全相同概念。
# HTO detects multiplexing collisions; scDblFinder uses transcriptomes.
#
# 【Output / 输出】
#   demultiplex_summary.csv
#   demultiplex_sample_summary.csv
#   demultiplex_applicability.csv
#   demultiplex_classification_counts.csv（仅HTO模式）
#   demultiplex_classification_overview.png
#   M05_demux_list.rds
#
# 【Check after running / 运行后检查】
# 检查每个HTO标签、Singlet比例和标签的真实生物学含义。
# Inspect HTO identities and singlet/doublet proportions.
# ============================================================

if (RUN_M05_DEMULTIPLEX) {

  MODULE <- "M05_demultiplex"
  MDIR <- module_dir(MODULE)

  seurat_list <- load_checkpoint(
    "M04_seurat_list"
  )

  # 保存HTODemux过滤前的真实classification计数。若配置为NONE或外部
  # metadata模式，则保留明确的适用性记录，而不是伪造HTO结果。
  # Preserve real HTODemux classifications before optional singlet filtering.
  # Non-HTO modes receive an explicit applicability record instead of a
  # simulated HTO chart.
  demux_classification_rows <- list()


  if (MULTIPLEXING_MODE == "NONE") {

    message(
      "MULTIPLEXING_MODE='NONE'：跳过HTO解复用。"
    )

    demux_list <- seurat_list


  } else if (
    MULTIPLEXING_MODE == "HTO"
  ) {

    demux_list <- list()

    for (
      sample_id in names(
        seurat_list
      )
    ) {

      obj <- seurat_list[[sample_id]]

      if (!"HTO" %in% SeuratObject::Assays(obj)) {

        safe_stop_module(
          MODULE,
          paste0(
            sample_id,
            " 没有HTO assay。\n",
            "如果这是多组织/多样本pool，不能直接进入组织筛选。\n",
            "请检查feature-barcode数据是否另存文件，或使用EXTERNAL_CELL_METADATA。"
          )
        )
      }


      DefaultAssay(obj) <- "HTO"

      obj <- NormalizeData(
        obj,
        normalization.method = "CLR",
        margin = 2,
        verbose = FALSE
      )

      obj <- HTODemux(
        obj,
        assay = "HTO",
        positive.quantile =
          HTO_POSITIVE_QUANTILE
      )


      classification_values <- as.character(
        obj$HTO_classification.global
      )

      classification_values[
        is.na(classification_values) |
          classification_values == ""
      ] <- "Unclassified"

      demux_classification_rows[[
        length(demux_classification_rows) + 1L
      ]] <- data.frame(
        sample_id = sample_id,
        classification = classification_values,
        stringsAsFactors = FALSE
      ) |>
        dplyr::count(
          sample_id,
          classification,
          name = "n_cells"
        )


      # HTO_maxID通常是最高的hashtag标签。
      raw_label <- as.character(
        obj$HTO_maxID
      )

      mapped_label <- raw_label

      if (length(HTO_LABEL_MAP) > 0L) {

        mapped <- unname(
          HTO_LABEL_MAP[
            raw_label
          ]
        )

        keep_original <- (
          is.na(mapped) |
          mapped == ""
        )

        mapped[
          keep_original
        ] <- raw_label[
          keep_original
        ]

        mapped_label <- mapped
      }


      if (HTO_LABEL_ROLE == "TISSUE") {

        # If labels are codes and no mapping was supplied, do not guess tissue identity.
        # If labels already contain the target tissue name, mapping is unnecessary.
        label_values <- unique(mapped_label[!is.na(mapped_label)])
        target_seen <- any(
          tolower(label_values) == tolower(TARGET_TISSUE)
        )

        if (
          length(HTO_LABEL_MAP) == 0L &&
          !target_seen
        ) {
          safe_stop_module(
            MODULE,
            paste0(
              sample_id,
              ": HTO_LABEL_ROLE='TISSUE', but TARGET_TISSUE='",
              TARGET_TISSUE,
              "' is not present among HTO labels and HTO_LABEL_MAP is empty.\n",
              "v6 will not guess tissue identities. Observed labels: ",
              paste(label_values, collapse = ", ")
            )
          )
        }

        obj[[TISSUE_COL]] <- mapped_label
      } else if (
        HTO_LABEL_ROLE == "SAMPLE"
      ) {
        obj$demux_sample <- mapped_label
      } else if (
        HTO_LABEL_ROLE == "DONOR"
      ) {
        obj$donor <- mapped_label
      } else {
        obj$hto_label <- mapped_label
      }


      if (KEEP_HTO_SINGLETS_ONLY) {

        keep_cells <- rownames(
          obj[[]]
        )[
          obj$HTO_classification.global ==
            "Singlet"
        ]

        obj <- subset(
          obj,
          cells = keep_cells
        )
      }

      demux_list[[sample_id]] <- obj
    }


  } else if (
    MULTIPLEXING_MODE ==
      "EXTERNAL_CELL_METADATA"
  ) {

    if (
      !nzchar(
        EXTERNAL_CELL_METADATA_FILE
      ) ||
      !file.exists(
        EXTERNAL_CELL_METADATA_FILE
      )
    ) {
      safe_stop_module(
        MODULE,
        "请提供EXTERNAL_CELL_METADATA_FILE。"
      )
    }

    cell_meta <- data.table::fread(
      EXTERNAL_CELL_METADATA_FILE,
      data.table = FALSE
    )

    demux_list <- seurat_list

    # 外部metadata具体列名可能因实验而异。
    # 为避免错误自动匹配，这里只检查常用字段。
    required <- c(
      SAMPLE_ID_COL,
      TISSUE_COL
    )

    missing <- setdiff(
      required,
      colnames(cell_meta)
    )

    if (length(missing) > 0L) {
      safe_stop_module(
        MODULE,
        paste0(
          "外部cell metadata缺少：",
          paste(
            missing,
            collapse = ", "
          )
        )
      )
    }

    # 真正barcode字段允许cell或cell_barcode。
    barcode_col <- if (
      "cell" %in%
        colnames(cell_meta)
    ) {
      "cell"
    } else if (
      "cell_barcode" %in%
        colnames(cell_meta)
    ) {
      "cell_barcode"
    } else {
      NA_character_
    }

    if (is.na(barcode_col)) {
      safe_stop_module(
        MODULE,
        "外部cell metadata必须包含cell或cell_barcode列。"
      )
    }


    for (
      sample_id in names(
        demux_list
      )
    ) {

      obj <- demux_list[[sample_id]]

      cm <- cell_meta[
        cell_meta[[SAMPLE_ID_COL]] ==
          sample_id,
        ,
        drop = FALSE
      ]

      idx <- match(
        colnames(obj),
        cm[[barcode_col]]
      )

      # 如果对象有sample前缀，而外部表只有原barcode，
      # 再尝试去除sample_id_。
      if (all(is.na(idx))) {

        stripped <- sub(
          paste0(
            "^",
            sample_id,
            "_"
          ),
          "",
          colnames(obj)
        )

        idx <- match(
          stripped,
          cm[[barcode_col]]
        )
      }

      if (anyNA(idx)) {
        safe_stop_module(
          MODULE,
          paste0(
            sample_id,
            ": ",
            sum(is.na(idx)),
            " cells are missing from external metadata / ",
            "个cell无法匹配外部metadata。"
          )
        )
      }

      obj[[TISSUE_COL]] <- cm[[TISSUE_COL]][
        idx
      ]

      demux_list[[sample_id]] <- obj
    }


  } else {

    safe_stop_module(
      MODULE,
      "不支持的MULTIPLEXING_MODE。"
    )
  }


  save_checkpoint(
    demux_list,
    "M05_demux_list"
  )


  demux_summary <- lapply(
    names(demux_list),
    function(sample_id) {

      obj <- demux_list[[sample_id]]

      metadata <- obj[[]]

      data.frame(
        sample_id = sample_id,
        n_cells = ncol(obj),
        has_tissue = TISSUE_COL %in%
          colnames(metadata),
        stringsAsFactors = FALSE
      )
    }
  ) |>
    data.table::rbindlist(
      fill = TRUE
    )

  data.table::fwrite(
    demux_summary,
    file.path(
      MDIR,
      "tables",
      "demultiplex_summary.csv"
    )
  )


  # 一个prepared RDS可能包含多个真实样本，因此额外按cell-level sample_id
  # 汇总。该表用于NONE/外部metadata模式的可视化，避免把RDS容器误画成
  # 一个生物学样本。
  # A prepared RDS may contain several biological samples. Summarizing the
  # cell-level sample field prevents an input container from being displayed
  # as though it were one biological replicate.
  demux_sample_summary_rows <- lapply(
    names(demux_list),
    function(object_id) {

      obj <- demux_list[[object_id]]
      metadata <- obj[[]]

      if (
        SAMPLE_ID_COL %in%
          colnames(metadata)
      ) {

        metadata |>
          tibble::rownames_to_column(
            "cell"
          ) |>
          dplyr::transmute(
            sample_id = as.character(
              .data[[SAMPLE_ID_COL]]
            ),
            condition = if (
              CONDITION_COL %in%
                colnames(metadata)
            ) {
              as.character(
                .data[[CONDITION_COL]]
              )
            } else {
              "Unspecified"
            }
          ) |>
          dplyr::count(
            sample_id,
            condition,
            name = "n_cells"
          )

      } else {

        data.frame(
          sample_id = object_id,
          condition = "Unspecified",
          n_cells = ncol(obj),
          stringsAsFactors = FALSE
        )
      }
    }
  )

  demux_sample_summary <- data.table::rbindlist(
    demux_sample_summary_rows,
    fill = TRUE
  ) |>
    as.data.frame() |>
    dplyr::group_by(
      sample_id,
      condition
    ) |>
    dplyr::summarise(
      n_cells = sum(
        .data$n_cells
      ),
      .groups = "drop"
    )

  data.table::fwrite(
    demux_sample_summary,
    file.path(
      MDIR,
      "tables",
      "demultiplex_sample_summary.csv"
    )
  )


  demultiplex_applicability <- data.frame(
    multiplexing_mode = MULTIPLEXING_MODE,
    hto_demultiplexing_applicable =
      MULTIPLEXING_MODE == "HTO",
    interpretation = if (
      MULTIPLEXING_MODE == "HTO"
    ) {
      "HTO classifications were calculated from the HTO assay."
    } else if (
      MULTIPLEXING_MODE == "EXTERNAL_CELL_METADATA"
    ) {
      "Cell identities were imported from externally validated metadata; HTODemux was not run."
    } else {
      "Each library already represents a known specimen; HTO demultiplexing is not applicable."
    },
    stringsAsFactors = FALSE
  )

  data.table::fwrite(
    demultiplex_applicability,
    file.path(
      MDIR,
      "tables",
      "demultiplex_applicability.csv"
    )
  )


  if (length(demux_classification_rows) > 0L) {

    demux_classification <- data.table::rbindlist(
      demux_classification_rows,
      fill = TRUE
    ) |>
      as.data.frame()

    data.table::fwrite(
      demux_classification,
      file.path(
        MDIR,
        "tables",
        "demultiplex_classification_counts.csv"
      )
    )

    p_demultiplex <- ggplot2::ggplot(
      demux_classification,
      ggplot2::aes(
        x = .data$sample_id,
        y = .data$n_cells,
        fill = .data$classification
      )
    ) +
      ggplot2::geom_col(
        position = "fill",
        width = 0.72,
        colour = "white",
        linewidth = 0.25
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Singlet" = "#54A76B",
          "Doublet" = "#F0A23B",
          "Negative" = "#9E9E9E",
          "Unclassified" = "#6C75A8"
        ),
        na.value = "#6C75A8",
        name = "HTO classification"
      ) +
      ggplot2::scale_y_continuous(
        labels = function(x) {
          paste0(
            round(
              100 * x
            ),
            "%"
          )
        },
        expand = ggplot2::expansion(
          mult = c(
            0,
            0.04
          )
        )
      ) +
      ggplot2::labs(
        title = "HTO demultiplexing classification",
        subtitle = "Proportions are calculated before optional singlet-only filtering",
        x = NULL,
        y = "Cell proportion",
        caption = "HTO doublets describe hashtag collisions and are not identical to RNA-based doublets."
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

  } else {

    # 无HTO时仍输出一个数据驱动的模块图：柱高为M05保留的真实细胞数，
    # subtitle明确说明为什么没有Singlet/Doublet/Negative分类。
    # When HTO is not applicable, the bars report actual retained cell counts
    # and the subtitle records why no HTO classes exist.
    p_demultiplex <- ggplot2::ggplot(
      demux_sample_summary,
      ggplot2::aes(
        x = stats::reorder(
          .data$sample_id,
          .data$n_cells
        ),
        y = .data$n_cells,
        fill = .data$condition
      )
    ) +
      ggplot2::geom_col(
        width = 0.68
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = format(
            .data$n_cells,
            big.mark = ",",
            scientific = FALSE,
            trim = TRUE
          )
        ),
        hjust = -0.12,
        size = 3.1
      ) +
      ggplot2::coord_flip(
        clip = "off"
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(
          mult = c(
            0,
            0.18
          )
        )
      ) +
      ggplot2::labs(
        title = "Cell identity assignment audit",
        subtitle = demultiplex_applicability$interpretation,
        x = NULL,
        y = "Cells retained in M05",
        fill = "Condition",
        caption = paste0(
          "Configured multiplexing mode: ",
          MULTIPLEXING_MODE
        )
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )
  }

  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "demultiplex_classification_overview.png"
    ),
    p_demultiplex,
    width = 8.5,
    height = max(
      5.5,
      0.38 *
        if (
          length(demux_classification_rows) > 0L
        ) {
          dplyr::n_distinct(
            demux_classification$sample_id
          )
        } else {
          nrow(demux_sample_summary)
        } +
        3.2
    ),
    dpi = 300
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M06. TARGET-TISSUE SUBSETTING AND MERGING / 目标组织筛选与合并
#
# 【Module scope / 模块范围】
# Tissue identity is assigned from the configured evidence source, optionally used
# to retain TARGET_TISSUE, and sample objects are then merged without changing raw
# counts. The original sample and replicate fields remain attached to every cell.
# 本模块依据配置的证据来源确定组织身份，可选择仅保留TARGET_TISSUE，随后在不改变
# 原始计数的情况下合并样本对象。每个细胞仍保留原始sample和replicate字段。
#
# 【Quality gate / 质量门槛】
# Cell counts must be reviewed before and after subsetting for every replicate.
# Loss of an entire replicate or a strongly asymmetric reduction may indicate an
# incorrect tissue map and must be resolved before QC or group comparison.
# 必须逐个生物学重复比较组织筛选前后的细胞数。若某个重复完全丢失或不同重复的
# 删除比例严重不对称，可能提示组织映射错误，应在QC和组间比较前解决。
#
# 【Purpose / 目的】
# 如果输入包含多个组织，只保留TARGET_TISSUE；然后合并各样本。
# Keep the target tissue when necessary, then merge samples.
#
# 【两种常见情况 / Two common cases】
# A. 自己实验每个sample本身就是Liver：
#    Each library is already liver-only:
#       TISSUE_ASSIGNMENT_MODE = "SAMPLE_METADATA"
#
# B. 一个library混有多种组织：
#    One library contains multiple tissues:
#       TISSUE_ASSIGNMENT_MODE = "HTO" or "EXTERNAL_CELL_METADATA"
#
# 【Output / 输出】
#   cells_after_tissue_subset.csv
#   tissue_subset_retention.csv
#   tissue_subset_cell_retention.png
#   M06_target_tissue_raw.rds
#
# 【Check after running / 运行后检查】
# 每个biological replicate筛选后是否仍保留足够细胞？
# Does every biological replicate retain a reasonable number of cells?
# ============================================================

if (RUN_M06_TISSUE_SUBSET) {

  MODULE <- "M06_tissue_subset"
  MDIR <- module_dir(MODULE)

  obj_list <- load_checkpoint(
    "M05_demux_list"
  )


  # 在任何组织筛选发生之前按真实实验设计字段统计细胞数。
  # Count cells by the configured biological-design fields before any tissue
  # subsetting so that sample loss and asymmetric retention remain auditable.
  cells_before_subset <- lapply(
    obj_list,
    function(obj) {

      metadata <- obj[[]]

      required_columns <- c(
        SAMPLE_ID_COL,
        BIOLOGICAL_REPLICATE_COL,
        CONDITION_COL
      )

      if (
        !all(
          required_columns %in%
            colnames(metadata)
        )
      ) {
        return(NULL)
      }

      metadata |>
        tibble::rownames_to_column(
          "cell"
        ) |>
        dplyr::transmute(
          sample_id = as.character(
            .data[[SAMPLE_ID_COL]]
          ),
          biological_replicate = as.character(
            .data[[BIOLOGICAL_REPLICATE_COL]]
          ),
          condition = as.character(
            .data[[CONDITION_COL]]
          )
        ) |>
        dplyr::count(
          sample_id,
          biological_replicate,
          condition,
          name = "n_cells_before"
        )
    }
  )

  cells_before_subset <- Filter(
    Negate(is.null),
    cells_before_subset
  )

  if (length(cells_before_subset) == 0L) {
    safe_stop_module(
      MODULE,
      paste0(
        "Cannot summarize pre-subset retention because required metadata are missing: ",
        paste(
          c(
            SAMPLE_ID_COL,
            BIOLOGICAL_REPLICATE_COL,
            CONDITION_COL
          ),
          collapse = ", "
        ),
        " / 组织筛选前对象缺少关键metadata，无法计算保留率。"
      )
    )
  }

  cells_before_subset <- data.table::rbindlist(
    cells_before_subset,
    fill = TRUE
  ) |>
    as.data.frame() |>
    dplyr::group_by(
      sample_id,
      biological_replicate,
      condition
    ) |>
    dplyr::summarise(
      n_cells_before = sum(
        .data$n_cells_before
      ),
      .groups = "drop"
    )


  filtered_list <- list()

  for (
    sample_id in names(
      obj_list
    )
  ) {

    obj <- obj_list[[sample_id]]


    if (FILTER_TARGET_TISSUE) {

      if (
        TISSUE_ASSIGNMENT_MODE ==
          "SAMPLE_METADATA"
      ) {

        tissue_value <- unique(
          as.character(
            obj[[]][[TISSUE_COL]]
          )
        )

        tissue_value <- tissue_value[
          !is.na(tissue_value)
        ]

        if (
          length(tissue_value) == 1L &&
          identical(
            tolower(
              tissue_value
            ),
            tolower(
              TARGET_TISSUE
            )
          )
        ) {

          # 整个sample都是目标组织，不需要筛细胞。
          obj[[TISSUE_COL]] <- tissue_value

        } else {

          # sample本身不是目标组织。
          next
        }


      } else if (
        TISSUE_ASSIGNMENT_MODE %in%
          c(
            "HTO",
            "EXTERNAL_CELL_METADATA"
          )
      ) {

        if (!TISSUE_COL %in%
            colnames(obj[[]])) {
          safe_stop_module(
            MODULE,
            paste0(
              sample_id,
              " 缺少cell-level tissue标签。"
            )
          )
        }

        keep <- colnames(obj)[
          tolower(
            as.character(
              obj[[]][[TISSUE_COL]]
            )
          ) ==
            tolower(
              TARGET_TISSUE
            )
        ]

        if (length(keep) == 0L) {

          warning(
            sample_id,
            " 没有识别到目标组织 ",
            TARGET_TISSUE
          )

          next
        }

        obj <- subset(
          obj,
          cells = keep
        )


      } else if (
        TISSUE_ASSIGNMENT_MODE ==
        "NONE"
      ) {
        safe_stop_module(
          MODULE,
          paste0(
            "FILTER_TARGET_TISSUE=TRUE requires a tissue-assignment mode. ",
            "要求筛选目标组织时，TISSUE_ASSIGNMENT_MODE不能为'NONE'。"
          )
        )
      } else {

        safe_stop_module(
          MODULE,
          "不支持的TISSUE_ASSIGNMENT_MODE。"
        )
      }
    }


    filtered_list[[sample_id]] <- obj
  }


  if (length(filtered_list) == 0L) {
    safe_stop_module(
      MODULE,
      "组织筛选后没有保留任何sample/cell。"
    )
  }


  if (length(filtered_list) == 1L) {

    seurat_raw <- filtered_list[[1L]]

  } else {

    seurat_raw <- Reduce(
      f = function(x, y) {
        merge(
          x,
          y,
          merge.data = FALSE
        )
      },
      x = filtered_list
    )
  }


  # 确保关键metadata存在。
  required_cell_metadata <- c(
    SAMPLE_ID_COL,
    BIOLOGICAL_REPLICATE_COL,
    CONDITION_COL
  )

  missing_cell_metadata <- setdiff(
    required_cell_metadata,
    colnames(
      seurat_raw[[]]
    )
  )

  if (
    length(
      missing_cell_metadata
    ) > 0L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "合并对象缺少关键metadata：",
        paste(
          missing_cell_metadata,
          collapse = ", "
        )
      )
    )
  }


  summary_table <- seurat_raw[[]] |>
    tibble::rownames_to_column(
      "cell"
    ) |>
    dplyr::count(
      .data[[SAMPLE_ID_COL]],
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]],
      name = "n_cells"
    )

  data.table::fwrite(
    summary_table,
    file.path(
      MDIR,
      "tables",
      "cells_after_tissue_subset.csv"
    )
  )


  cells_after_subset <- summary_table |>
    dplyr::transmute(
      sample_id = as.character(
        .data[[SAMPLE_ID_COL]]
      ),
      biological_replicate = as.character(
        .data[[BIOLOGICAL_REPLICATE_COL]]
      ),
      condition = as.character(
        .data[[CONDITION_COL]]
      ),
      n_cells_after = as.numeric(
        .data$n_cells
      )
    ) |>
    dplyr::group_by(
      sample_id,
      biological_replicate,
      condition
    ) |>
    dplyr::summarise(
      n_cells_after = sum(
        .data$n_cells_after
      ),
      .groups = "drop"
    )

  tissue_subset_retention <- dplyr::full_join(
    cells_before_subset,
    cells_after_subset,
    by = c(
      "sample_id",
      "biological_replicate",
      "condition"
    )
  ) |>
    dplyr::mutate(
      n_cells_before = dplyr::coalesce(
        .data$n_cells_before,
        0
      ),
      n_cells_after = dplyr::coalesce(
        .data$n_cells_after,
        0
      ),
      retention_percent = dplyr::if_else(
        .data$n_cells_before > 0,
        100 *
          .data$n_cells_after /
          .data$n_cells_before,
        NA_real_
      )
    ) |>
    dplyr::arrange(
      .data$condition,
      .data$sample_id,
      .data$biological_replicate
    )

  data.table::fwrite(
    tissue_subset_retention,
    file.path(
      MDIR,
      "tables",
      "tissue_subset_retention.csv"
    )
  )


  retention_plot_data <- tissue_subset_retention |>
    dplyr::mutate(
      sample_label = paste0(
        .data$sample_id,
        " [",
        .data$biological_replicate,
        "]"
      )
    ) |>
    dplyr::select(
      sample_label,
      condition,
      n_cells_before,
      n_cells_after
    ) |>
    tidyr::pivot_longer(
      cols = c(
        n_cells_before,
        n_cells_after
      ),
      names_to = "stage",
      values_to = "n_cells"
    ) |>
    dplyr::mutate(
      stage = factor(
        .data$stage,
        levels = c(
          "n_cells_before",
          "n_cells_after"
        ),
        labels = c(
          "Before tissue subset",
          "After tissue subset"
        )
      )
    )

  p_tissue_retention <- ggplot2::ggplot(
    retention_plot_data,
    ggplot2::aes(
      x = .data$sample_label,
      y = .data$n_cells,
      fill = .data$stage
    )
  ) +
    ggplot2::geom_col(
      position = ggplot2::position_dodge(
        width = 0.76
      ),
      width = 0.68
    ) +
    ggplot2::facet_wrap(
      ggplot2::vars(
        condition
      ),
      scales = "free_x"
    ) +
    ggplot2::scale_fill_manual(
      values = c(
        "Before tissue subset" = "#9ECAE1",
        "After tissue subset" = "#2B6DA8"
      ),
      name = NULL
    ) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.08
        )
      )
    ) +
    ggplot2::labs(
      title = "Cell retention after target-tissue selection",
      subtitle = if (
        FILTER_TARGET_TISSUE
      ) {
        paste0(
          "Target tissue: ",
          TARGET_TISSUE,
          "; assignment mode: ",
          TISSUE_ASSIGNMENT_MODE
        )
      } else {
        "No tissue filtering was requested; before and after counts should be identical"
      },
      x = NULL,
      y = "Number of cells",
      caption = "A missing replicate or strongly asymmetric retention should be resolved before QC."
    ) +
    ggplot2::theme_minimal(
      base_size = 11
    ) +
    ggplot2::theme(
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1
      ),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )

  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "tissue_subset_cell_retention.png"
    ),
    p_tissue_retention,
    width = max(
      9,
      0.58 *
        dplyr::n_distinct(
          retention_plot_data$sample_label
        ) +
        4
    ),
    height = 6.5,
    dpi = 300
  )

  print(summary_table)

  save_checkpoint(
    seurat_raw,
    "M06_target_tissue_raw"
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M07. QUALITY CONTROL / 质量控制
#
# 【Module scope / 模块范围】
# Per-cell RNA quality metrics are calculated, fixed and/or adaptive failure flags
# are recorded, and filtering is applied only when APPLY_QC_FILTER = TRUE. Both
# pre-filter and post-filter summaries are retained for traceability.
# 本模块计算逐细胞RNA质量指标，记录固定阈值和/或自适应异常值标记，并仅在
# APPLY_QC_FILTER = TRUE时执行过滤。过滤前后汇总均会保留，以便追溯。
#
# 【Interpretation boundary / 解读边界】
# QC removes cells with profiles consistent with technical failure; it is not a
# method for selecting desired biological phenotypes. Thresholds should be judged
# from per-sample distributions, retained cell numbers and tissue-specific biology,
# especially for RNA-rich hepatocytes and stressed or injured populations.
# QC用于排除符合技术失败特征的细胞，不是选择“期望”生物学表型的方法。阈值应结合
# 各样本分布、保留细胞数和组织特异性生物学判断，尤其需谨慎对待RNA含量较高的
# hepatocyte以及应激或损伤相关细胞群。
#
# 【OSCA对应 / OSCA connection】
# OSCA basic - Quality Control.
# 主要关注三个常见指标：
# Three common metrics are emphasized:
#   library size      ≈ nCount_RNA
#   detected features ≈ nFeature_RNA
#   mitochondrial %   ≈ percent.mt
#
# 【本脚本为什么默认OSCA模式 / Why OSCA mode is default】
# OSCA推荐通过样本分布识别outlier，而不是跨所有数据机械套统一阈值。
# OSCA emphasizes distribution-aware outlier detection rather than
# blindly applying one threshold to every experiment.
#
# computeRnaQcMetrics():
#   计算每个细胞的总counts、检测基因数和线粒体比例。
#   Calculates library size, detected features and mitochondrial fraction.
#
# suggestRnaQcThresholds():
#   根据QC指标分布建议自适应异常值阈值。
#   Suggests adaptive QC thresholds from the observed QC distributions.
#
# batch = sample_id:
#   在每个样本内部判断异常值，降低测序深度差异造成的误删风险。
#   Detects outliers within sample to reduce inappropriate filtering
#   caused by sample-specific sequencing depth.
#
# 【三种模式 / Three modes】
#   OSCA   = adaptive QC only
#   FIXED  = user-defined hard thresholds only
#   HYBRID = both rules must pass
#
# 【重要假设 / Important assumption】
# QC指标的极端值主要来自技术质量，而不是某个真实细胞类型本身。
# Extreme QC values are assumed to reflect technical quality rather
# than genuine biology. Always inspect diagnostic plots.
#
# 【Output / 输出】
#   OSCA_QC_flags.csv
#   QC_summary.csv
#   QC_violin_before.png
#   QC_violin_after.png
#   M07_after_QC.rds
#
# 【Check after running / 运行后检查】
# 1. 每个sample删除比例 / fraction removed per sample
# 2. 是否某个sample几乎被删光 / whether one sample is over-filtered
# 3. 肝细胞高RNA含量是否被误当异常 / possible loss of RNA-rich hepatocytes
# ============================================================

if (RUN_M07_QC) {

  MODULE <- "M07_QC"
  MDIR <- module_dir(MODULE)

  seurat_raw <- load_checkpoint(
    "M06_target_tissue_raw"
  )

  DefaultAssay(
    seurat_raw
  ) <- "RNA"


  mitochondrial_genes <- grep(
    species_cfg$mito_pattern,
    rownames(seurat_raw),
    value = TRUE
  )

  if (
    length(
      mitochondrial_genes
    ) > 0L
  ) {

    seurat_raw[["percent.mt"]] <- PercentageFeatureSet(
      seurat_raw,
      features =
        mitochondrial_genes
    )

  } else {

    warning(
      "未识别到线粒体基因。请检查基因命名。"
    )

    seurat_raw$percent.mt <- NA_real_
  }


  ribosomal_genes <- grep(
    species_cfg$ribo_pattern,
    rownames(seurat_raw),
    value = TRUE
  )

  if (
    length(
      ribosomal_genes
    ) > 0L
  ) {

    seurat_raw[["percent.ribo"]] <- PercentageFeatureSet(
      seurat_raw,
      features =
        ribosomal_genes
    )

  } else {
    seurat_raw$percent.ribo <- NA_real_
  }


  # ---------- 固定阈值 ----------
  fixed_pass <- (
    seurat_raw$nFeature_RNA >=
      QC_MIN_FEATURES &
    seurat_raw$nFeature_RNA <=
      QC_MAX_FEATURES &
    seurat_raw$nCount_RNA >=
      QC_MIN_COUNTS &
    seurat_raw$nCount_RNA <=
      QC_MAX_COUNTS
  )

  if (
    any(
      is.finite(
        seurat_raw$percent.mt
      )
    )
  ) {
    fixed_pass <- (
      fixed_pass &
      seurat_raw$percent.mt <=
        QC_MAX_MT
    )
  }


  # ---------- OSCA QC ----------
  #
  # Current Bioconductor implementation using scrapper.
  # 使用当前Bioconductor scrapper API进行OSCA风格QC。
  #
  # The general strategy follows the current OSCA/scrapper
  # recommendation:
  #
  #   raw RNA counts
  #       ↓
  #   per-cell QC metrics
  #       ↓
  #   sample-aware adaptive thresholds
  #       ↓
  #   high-quality / low-quality classification
  #
  # 当前流程遵循OSCA/scrapper推荐思路：
  #
  #   原始RNA counts
  #       ↓
  #   计算每个cell的QC指标
  #       ↓
  #   在sample内部确定自适应阈值
  #       ↓
  #   判断cell是否通过QC
  #
  # IMPORTANT:
  #
  # computeRnaQcMetrics() expects a matrix-like count object,
  # not a SingleCellExperiment object.
  #
  # computeRnaQcMetrics()需要输入RNA count矩阵，
  # 不能直接输入SingleCellExperiment对象。
  # ----------------------------------------------------------


  # ----------------------------------------------------------
  # Join Seurat v5 RNA count layers if necessary
  # 如有必要，合并Seurat v5 RNA count layers
  # ----------------------------------------------------------
  #
  # A merged Seurat v5 object may contain:
  #
  # counts.Sample1
  # counts.Sample2
  # counts.Sample3
  #
  # rather than one single counts layer.
  #
  # 合并多个sample以后，
  # Seurat v5可能把不同sample保存在不同counts layer中。
  #
  # The helper function joins these layers before QC.
  #
  # 这里先合并这些layers，
  # 使后面的Bioconductor QC能够直接读取完整RNA count矩阵。
  #
  # This is NOT batch correction or integration.
  #
  # 注意：
  # JoinLayers不是批次校正，也不是数据整合，
  # 只是把同一个RNA assay中的count layers合并起来。
  # ----------------------------------------------------------

  seurat_raw <- prepare_for_sce(
    seurat_raw,
    assay = "RNA"
  )


  # ----------------------------------------------------------
  # Extract the raw RNA count matrix
  # 提取原始RNA count矩阵
  # ----------------------------------------------------------
  #
  # scrapper::computeRnaQcMetrics() expects:
  #
  # rows    = genes
  # columns = cells
  # values  = raw RNA counts
  #
  # 因此这里直接从Seurat RNA assay提取counts，
  # 而不是先转换成SingleCellExperiment。
  # ----------------------------------------------------------

  rna_counts_qc <- SeuratObject::LayerData(
    object = seurat_raw,
    assay = "RNA",
    layer = "counts"
  )


  # ----------------------------------------------------------
  # Safety check for RNA count matrix
  # 检查RNA count矩阵
  # ----------------------------------------------------------

  if (ncol(rna_counts_qc) != ncol(seurat_raw)) {

    safe_stop_module(
      MODULE,
      paste0(
        "RNA count matrix cell number does not match Seurat object. ",
        "RNA count矩阵cell数量与Seurat对象不一致。 ",
        "counts = ",
        ncol(rna_counts_qc),
        "; Seurat = ",
        ncol(seurat_raw)
      )
    )
  }


  # ----------------------------------------------------------
  # Define mitochondrial genes
  # 定义线粒体基因
  # ----------------------------------------------------------
  #
  # mitochondrial_genes was identified earlier according
  # to the species-specific mitochondrial gene pattern.
  #
  # mitochondrial_genes已经在前面根据物种对应的
  # 线粒体基因命名规则识别。
  #
  # scrapper accepts a logical vector identifying the
  # mitochondrial genes.
  #
  # scrapper可以直接使用logical vector标记
  # 哪些基因属于线粒体基因。
  # ----------------------------------------------------------

  is_mito <- rownames(
    rna_counts_qc
  ) %in% mitochondrial_genes


  # ----------------------------------------------------------
  # Calculate per-cell RNA QC metrics
  # 计算每个cell的RNA QC指标
  # ----------------------------------------------------------
  #
  # computeRnaQcMetrics() calculates metrics including:
  #
  # sum
  #   total RNA counts / library size
  #   总RNA counts
  #
  # detected
  #   number of detected genes
  #   检测到的基因数量
  #
  # subsets$Mito
  #   mitochondrial proportion
  #   线粒体reads/count比例
  #
  # These correspond to the core QC concepts used in OSCA.
  #
  # 这些正是OSCA QC中最核心的指标。
  # ----------------------------------------------------------

  qc_stats <- scrapper::computeRnaQcMetrics(
    rna_counts_qc,
    subsets = list(
      Mito = is_mito
    )
  )


  # ----------------------------------------------------------
  # Extract sample/batch information
  # 提取sample/batch信息
  # ----------------------------------------------------------
  #
  # QC_BATCH_COL is defined in Section 0.
  #
  # 默认：
  #
  # QC_BATCH_COL <- SAMPLE_ID_COL
  #
  # meaning that adaptive QC thresholds are estimated
  # separately for each biological library/sample.
  #
  # 即默认在每个sample内部确定自适应QC阈值。
  #
  # This is useful because sequencing depth and capture
  # efficiency may differ between samples.
  #
  # 不同sample的测序深度和捕获效率可能不同，
  # 因此sample-aware QC通常比所有cells使用同一个
  # adaptive threshold更加合理。
  # ----------------------------------------------------------

  cell_metadata <- seurat_raw[[]]


  if (!(QC_BATCH_COL %in% colnames(cell_metadata))) {

    safe_stop_module(
      MODULE,
      paste0(
        "QC batch column not found / ",
        "没有找到QC分组列: ",
        QC_BATCH_COL
      )
    )
  }


  batch_values <- as.character(
    cell_metadata[[QC_BATCH_COL]]
  )


  # ----------------------------------------------------------
  # Safety check
  # 安全检查
  # ----------------------------------------------------------

  if (length(batch_values) != ncol(seurat_raw)) {

    safe_stop_module(
      MODULE,
      paste0(
        "QC batch vector length does not match cell number. ",
        "QC batch vector长度与cell数量不一致。 ",
        "batch length = ",
        length(batch_values),
        "; cells = ",
        ncol(seurat_raw)
      )
    )
  }


  if (anyNA(batch_values)) {

    safe_stop_module(
      MODULE,
      paste0(
        "QC batch column contains NA values / ",
        "QC分组列存在NA: ",
        QC_BATCH_COL
      )
    )
  }


  # Convert to factor for block-aware threshold estimation.
  # 转换成factor，用于block-aware QC。

  batch_values <- factor(
    batch_values
  )


  # ----------------------------------------------------------
  # Suggest adaptive QC thresholds
  # 计算自适应QC阈值
  # ----------------------------------------------------------
  #
  # Current scrapper uses the argument:
  #
  # block =
  #
  # rather than the older:
  #
  # batch =
  #
  # 当前scrapper API使用block参数。
  #
  # A separate threshold is estimated for each sample/block.
  #
  # 每个sample分别计算阈值。
  #
  # By default, scrapper uses MAD-based adaptive thresholds.
  #
  # 默认使用基于MAD的自适应异常值判断，
  # 与OSCA推荐的QC思想一致。
  # ----------------------------------------------------------

  qc_thresholds <- scrapper::suggestRnaQcThresholds(
    qc_stats,
    block = batch_values
  )


  # ----------------------------------------------------------
  # Apply adaptive QC thresholds
  # 应用自适应QC阈值
  # ----------------------------------------------------------
  #
  # filterRnaQcMetrics() returns:
  #
  # TRUE  = cell passes adaptive QC
  # FALSE = cell fails adaptive QC
  #
  # filterRnaQcMetrics()直接返回每个cell是否通过QC。
  # ----------------------------------------------------------

  osca_pass <- scrapper::filterRnaQcMetrics(
    qc_thresholds,
    qc_stats,
    block = batch_values
  )


  # ----------------------------------------------------------
  # Safety check for final adaptive QC vector
  # 检查最终OSCA QC结果
  # ----------------------------------------------------------

  if (length(osca_pass) != ncol(seurat_raw)) {

    safe_stop_module(
      MODULE,
      paste0(
        "OSCA QC result length does not match cell number. ",
        "OSCA QC结果长度与cell数量不一致。"
      )
    )
  }



  # ---------- 最终规则 ----------
  if (QC_MODE == "OSCA") {
    qc_pass <- osca_pass
  } else if (
    QC_MODE == "FIXED"
  ) {
    qc_pass <- fixed_pass
  } else if (
    QC_MODE == "HYBRID"
  ) {
    qc_pass <- (
      osca_pass &
      fixed_pass
    )
  } else {
    safe_stop_module(
      MODULE,
      "QC_MODE必须是OSCA/FIXED/HYBRID。"
    )
  }


  seurat_raw$qc_pass <- qc_pass
  seurat_raw$qc_pass_osca <- osca_pass
  seurat_raw$qc_pass_fixed <- fixed_pass


  # ============================================================
  # SAVE OSCA QC METRICS AND FLAGS
  # 保存OSCA QC指标和筛选结果
  # ============================================================
  #
  # In the current scrapper workflow:
  #
  #   qc_stats
  #     contains per-cell QC metrics
  #     保存每个cell的QC指标
  #
  #   qc_thresholds
  #     contains the adaptive thresholds
  #     保存自适应QC阈值
  #
  #   osca_pass
  #     TRUE  = pass adaptive QC
  #     FALSE = fail adaptive QC
  #
  # The older scuttle workflow generated an object called
  # qc_reasons. That object no longer exists in the current
  # scrapper implementation.
  #
  # 旧版scuttle流程会生成qc_reasons，
  # 但当前scrapper流程已经不再产生这个对象。
  #
  # Therefore we explicitly combine qc_stats and osca_pass
  # into a cell-level QC table.
  #
  # 因此这里把qc_stats与osca_pass合并，
  # 生成新的cell-level QC结果表。
  # ============================================================

  if (SAVE_QC_FLAGS) {

    qr <- as.data.frame(
      qc_stats
    )

    # ----------------------------------------------------------
    # Add cell identifiers
    # 添加cell ID
    # ----------------------------------------------------------

    qr$cell <- colnames(
      seurat_raw
    )


    # ----------------------------------------------------------
    # Add sample/batch information
    # 添加sample/batch信息
    # ----------------------------------------------------------

    qr$qc_block <- as.character(
      batch_values
    )


    # ----------------------------------------------------------
    # Add adaptive QC result
    # 添加OSCA自适应QC结果
    # ----------------------------------------------------------

    qr$qc_pass_osca <- as.logical(
      osca_pass
    )


    # ----------------------------------------------------------
    # Add fixed-threshold QC result
    # 添加固定阈值QC结果
    # ----------------------------------------------------------

    qr$qc_pass_fixed <- as.logical(
      fixed_pass
    )


    # ----------------------------------------------------------
    # Add final QC result
    # 添加最终QC结果
    # ----------------------------------------------------------
    #
    # qc_pass depends on QC_MODE:
    #
    # OSCA
    #   adaptive QC only
    #
    # FIXED
    #   fixed thresholds only
    #
    # HYBRID
    #   both must pass
    #
    # qc_pass由前面设置的QC_MODE决定。
    # ----------------------------------------------------------

    qr$qc_pass_final <- as.logical(
      qc_pass
    )


    # ----------------------------------------------------------
    # Safety check
    # 安全检查
    # ----------------------------------------------------------

    if (nrow(qr) != ncol(seurat_raw)) {

      safe_stop_module(
        MODULE,
        paste0(
          "QC table row number does not match cell number. ",
          "QC结果表行数与cell数量不一致。 ",
          "QC rows = ",
          nrow(qr),
          "; cells = ",
          ncol(seurat_raw)
        )
      )
    }


    # ----------------------------------------------------------
    # Save cell-level QC table
    # 保存cell-level QC结果
    # ----------------------------------------------------------

    data.table::fwrite(
      qr,
      file.path(
        MDIR,
        "tables",
        "OSCA_QC_flags.csv"
      )
    )
  }


  p_before <- VlnPlot(
    seurat_raw,
    features = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt"
    ),
    group.by = SAMPLE_ID_COL,
    pt.size = 0,
    ncol = 3
  ) &
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "QC_violin_before.png"
    ),
    p_before,
    width = 15,
    height = 6,
    dpi = 300
  )


  if (APPLY_QC_FILTER) {

    keep_cells <- colnames(
      seurat_raw
    )[
      seurat_raw$qc_pass
    ]

    seurat_qc <- subset(
      seurat_raw,
      cells = keep_cells
    )

  } else {

    seurat_qc <- seurat_raw
  }


  p_after <- VlnPlot(
    seurat_qc,
    features = c(
      "nCount_RNA",
      "nFeature_RNA",
      "percent.mt"
    ),
    group.by = SAMPLE_ID_COL,
    pt.size = 0,
    ncol = 3
  ) &
    theme(
      axis.text.x =
        element_text(
          angle = 45,
          hjust = 1
        )
    )

  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "QC_violin_after.png"
    ),
    p_after,
    width = 15,
    height = 6,
    dpi = 300
  )


  # ----------------------------------------------------------
  # Seurat v5 layer compatibility before checkpoint
  # 保存checkpoint前统一Seurat v5 RNA layers
  # ----------------------------------------------------------
  #
  # This is not batch correction or integration.
  # It only prevents downstream GetAssayData()/SCE conversion
  # failures caused by multiple Assay5 layers.
  #
  # 这不是批次校正或integration，只是避免多个Assay5 layers
  # 导致后续GetAssayData()/SingleCellExperiment转换失败。
  seurat_qc <- join_assay_layers_if_needed(
    seurat_qc,
    assay = "RNA",
    verbose = TRUE
  )

  save_checkpoint(
    seurat_qc,
    "M07_after_QC"
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M08. DOUBLET DETECTION / Doublet识别
#
# 【Module scope / 模块范围】
# Transcriptome-based doublet scores and classes are estimated within the declared
# sample structure. Detection can be run without removal, allowing sample-specific
# rates and marker patterns to be reviewed before cells are excluded.
# 本模块在已声明的样本结构内估计基于转录组的doublet分数和类别。可仅运行检测而不
# 删除，以便在排除细胞前检查各样本比例和marker模式。
#
# 【Interpretation boundary / 解读边界】
# A prediction is probabilistic and should be assessed alongside expected loading,
# library complexity and cluster markers. High RNA content alone is insufficient
# evidence, and unusually high sample-specific rates require technical review.
# doublet预测属于概率判定，应结合预期上样量、文库复杂度和cluster marker评估。
# 单纯RNA含量高不足以判定doublet；若某样本预测比例异常高，应回查实验和技术过程。
#
# 【概念 / Concept】
# Doublet表示一个液滴中捕获两个或更多细胞，导致混合表达谱。
# A doublet occurs when multiple cells are captured in one droplet,
# producing a mixed transcriptomic profile.
#
# 【为什么不能只看高UMI / Why not just use high UMI?】
# Hepatocytes等细胞本身RNA含量就高；高nCount_RNA不自动等于doublet。
# RNA-rich cells such as hepatocytes may naturally have high counts.
#
# scDblFinder:
#   利用真实细胞和模拟doublet的表达特征进行分类。
#   Uses transcriptomic information and simulated doublets for classification.
#
# samples = sample_id:
#   在不同文库内部估计，而不是把所有文库完全混在一起。
#   Performs sample-aware doublet detection across libraries.
#
# 【Output / 输出】
#   doublet_summary.csv
#   doublet_classification_by_sample.png
#   doublet_score_and_complexity.png
#   doublet_detection_status.png（仅在关闭scDblFinder时）
#   M08_singlets.rds
#
# 【Check after running / 运行后检查】
# 比较不同sample的doublet比例；异常高比例需要回看实验和QC。
# Compare doublet fractions across samples and investigate outliers.
# ============================================================

if (RUN_M08_DOUBLETS) {

  MODULE <- "M08_doublets"
  MDIR <- module_dir(MODULE)

  seurat_qc <- load_checkpoint(
    "M07_after_QC"
  )


  if (!RUN_SCDOUBLETFINDER) {

    seurat_singlet <- seurat_qc

    # 即使用户关闭doublet检测，也输出清晰的模块状态图，而不是让M08 figures
    # 文件夹保持空白。柱高来自当前对象的真实样本细胞数，不代表doublet结果。
    # When detection is disabled, report the actual cells entering M08 and state
    # explicitly that no singlet/doublet inference was performed.
    doublet_disabled_summary <- seurat_qc[[]] |>
      tibble::rownames_to_column(
        "cell"
      ) |>
      dplyr::transmute(
        sample_id = as.character(
          .data[[SAMPLE_ID_COL]]
        )
      ) |>
      dplyr::count(
        sample_id,
        name = "n_cells"
      )

    p_doublet_disabled <- ggplot2::ggplot(
      doublet_disabled_summary,
      ggplot2::aes(
        x = stats::reorder(
          .data$sample_id,
          .data$n_cells
        ),
        y = .data$n_cells
      )
    ) +
      ggplot2::geom_col(
        width = 0.68,
        fill = "#4C78A8"
      ) +
      ggplot2::coord_flip() +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(
          mult = c(
            0,
            0.08
          )
        )
      ) +
      ggplot2::labs(
        title = "Doublet-detection status",
        subtitle = "scDblFinder was disabled; no cells were classified or removed in M08",
        x = NULL,
        y = "Cells entering M08",
        caption = "Cell counts are descriptive and must not be interpreted as doublet evidence."
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.y = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        )
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "doublet_detection_status.png"
      ),
      p_doublet_disabled,
      width = 8.5,
      height = max(
        5.5,
        0.38 *
          nrow(
            doublet_disabled_summary
          ) +
          3.2
      ),
      dpi = 300
    )

  } else {

    set.seed(
      RANDOM_SEED
    )

    # Join possible Seurat v5 RNA layers before scDblFinder.
    # 在scDblFinder前先统一可能存在的RNA layers。
    seurat_qc <- prepare_for_sce(
      seurat_qc,
      assay = "RNA"
    )

    sce_doublet <- as.SingleCellExperiment(
      seurat_qc,
      assay = "RNA"
    )

    # scDblFinder expects the sample information to be available
    # in colData when samples is provided as a column name.
    SummarizedExperiment::colData(sce_doublet)[[SAMPLE_ID_COL]] <-
      as.character(
        seurat_qc[[SAMPLE_ID_COL]][, 1]
      )

    sce_doublet <- scDblFinder::scDblFinder(
      sce_doublet,
      samples = SAMPLE_ID_COL,
      BPPARAM =
        BiocParallel::SerialParam(
          RNGseed = RANDOM_SEED
        )
    )

    cd <- as.data.frame(
      SummarizedExperiment::colData(
        sce_doublet
      )
    )

    idx <- match(
      colnames(
        seurat_qc
      ),
      rownames(cd)
    )

    seurat_qc$scDblFinder.score <-
      cd$scDblFinder.score[
        idx
      ]

    seurat_qc$scDblFinder.class <-
      cd$scDblFinder.class[
        idx
      ]


    doublet_summary <- seurat_qc[[]] |>
      tibble::rownames_to_column(
        "cell"
      ) |>
      dplyr::count(
        .data[[SAMPLE_ID_COL]],
        scDblFinder.class,
        name = "n_cells"
      ) |>
      dplyr::group_by(
        .data[[SAMPLE_ID_COL]]
      ) |>
      dplyr::mutate(
        fraction =
          n_cells /
          sum(n_cells)
      ) |>
      dplyr::ungroup()

    data.table::fwrite(
      doublet_summary,
      file.path(
        MDIR,
        "tables",
        "doublet_summary.csv"
      )
    )


    # ========================================================
    # 1. SAMPLE-LEVEL CLASSIFICATION COMPOSITION
    # 1. 样本层面的doublet分类组成
    # ========================================================
    #
    # 该图使用所有细胞，并同时展示比例和doublet百分比。样本间异常高的
    # doublet比例可提示上样浓度、细胞聚集或文库质量需要回查，但不能单独
    # 证明算法预测必然正确。
    #
    # This plot uses every classified cell. An unusually high sample-specific
    # fraction is a technical-review flag, not proof that every prediction is
    # biologically correct.
    # ========================================================

    doublet_composition_plot <- doublet_summary |>
      dplyr::transmute(
        sample_id = as.character(
          .data[[SAMPLE_ID_COL]]
        ),
        classification = factor(
          as.character(
            .data$scDblFinder.class
          ),
          levels = c(
            "singlet",
            "doublet"
          ),
          labels = c(
            "Singlet",
            "Doublet"
          )
        ),
        n_cells = as.numeric(
          .data$n_cells
        ),
        fraction = as.numeric(
          .data$fraction
        )
      ) |>
      dplyr::mutate(
        sample_id = factor(
          .data$sample_id,
          levels = unique(
            .data$sample_id
          )
        )
      )

    p_doublet_composition <- ggplot2::ggplot(
      doublet_composition_plot,
      ggplot2::aes(
        x = .data$sample_id,
        y = .data$fraction,
        fill = .data$classification
      )
    ) +
      ggplot2::geom_col(
        width = 0.72,
        colour = "white",
        linewidth = 0.25
      ) +
      ggplot2::geom_text(
        data = doublet_composition_plot |>
          dplyr::filter(
            .data$classification ==
              "Doublet"
          ),
        ggplot2::aes(
          y = 1.035,
          label = paste0(
            formatC(
              100 * .data$fraction,
              format = "f",
              digits = 1
            ),
            "%"
          )
        ),
        colour = "#B3261E",
        fontface = "bold",
        size = 3.1,
        show.legend = FALSE
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Singlet" = "#4C78A8",
          "Doublet" = "#E45756"
        ),
        drop = FALSE,
        name = NULL
      ) +
      ggplot2::scale_y_continuous(
        labels = function(x) {
          paste0(
            round(
              100 * x
            ),
            "%"
          )
        },
        limits = c(
          0,
          1.08
        ),
        expand = ggplot2::expansion(
          mult = c(
            0,
            0
          )
        )
      ) +
      ggplot2::labs(
        title = "Predicted doublet composition by sample",
        subtitle = "Percentages above the red segments are scDblFinder-predicted doublet fractions",
        x = NULL,
        y = "Cell proportion",
        caption = "Predictions are sample-aware and should be reviewed with loading expectations and library complexity."
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "doublet_classification_by_sample.png"
      ),
      p_doublet_composition,
      width = max(
        8.5,
        0.65 *
          dplyr::n_distinct(
            doublet_composition_plot$sample_id
          ) +
          4
      ),
      height = 6.2,
      dpi = 300
    )


    # ========================================================
    # 2. SCORE AND LIBRARY-COMPLEXITY DIAGNOSTIC
    # 2. score与文库复杂度诊断
    # ========================================================
    #
    # 左图比较预测类别的score分布；右图检查预测doublet是否仅由高UMI或
    # 高feature数量驱动。右图为减少遮盖进行固定随机种子抽样，所有细胞的
    # 分类和上方比例图均保持完整。
    #
    # The score panel assesses class separation. The complexity panel tests
    # whether calls are concentrated among high-count/high-feature cells.
    # Deterministic subsampling affects only the scatter display.
    # ========================================================

    doublet_cell_diagnostics <- seurat_qc[[]] |>
      tibble::rownames_to_column(
        "cell"
      ) |>
      dplyr::transmute(
        cell = .data$cell,
        sample_id = as.character(
          .data[[SAMPLE_ID_COL]]
        ),
        classification = factor(
          as.character(
            .data$scDblFinder.class
          ),
          levels = c(
            "singlet",
            "doublet"
          ),
          labels = c(
            "Singlet",
            "Doublet"
          )
        ),
        score = as.numeric(
          .data$scDblFinder.score
        ),
        nCount_RNA = as.numeric(
          .data$nCount_RNA
        ),
        nFeature_RNA = as.numeric(
          .data$nFeature_RNA
        )
      ) |>
      dplyr::filter(
        !is.na(
          .data$classification
        ),
        is.finite(
          .data$score
        ),
        is.finite(
          .data$nCount_RNA
        ),
        is.finite(
          .data$nFeature_RNA
        )
      )

    set.seed(
      RANDOM_SEED
    )

    doublet_complexity_plot <- doublet_cell_diagnostics |>
      dplyr::filter(
        .data$nCount_RNA > 0,
        .data$nFeature_RNA > 0
      ) |>
      dplyr::group_by(
        sample_id,
        classification
      ) |>
      dplyr::group_modify(
        function(.x, .y) {

          n_keep <- min(
            nrow(.x),
            as.integer(
              DOUBLET_DIAGNOSTIC_MAX_CELLS_PER_SAMPLE_CLASS
            )
          )

          if (n_keep < nrow(.x)) {
            dplyr::slice_sample(
              .x,
              n = n_keep
            )
          } else {
            .x
          }
        }
      ) |>
      dplyr::ungroup()

    p_doublet_score <- ggplot2::ggplot(
      doublet_cell_diagnostics,
      ggplot2::aes(
        x = .data$classification,
        y = .data$score,
        fill = .data$classification
      )
    ) +
      ggplot2::geom_violin(
        scale = "width",
        trim = TRUE,
        alpha = 0.78,
        colour = NA
      ) +
      ggplot2::geom_boxplot(
        width = 0.16,
        outlier.shape = NA,
        fill = "white",
        colour = "#333333",
        linewidth = 0.35
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Singlet" = "#4C78A8",
          "Doublet" = "#E45756"
        ),
        guide = "none"
      ) +
      ggplot2::labs(
        title = "Prediction-score separation",
        x = NULL,
        y = "scDblFinder score"
      ) +
      ggplot2::theme_minimal(
        base_size = 10.5
      ) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        )
      )

    p_doublet_complexity <- ggplot2::ggplot(
      doublet_complexity_plot,
      ggplot2::aes(
        x = .data$nCount_RNA,
        y = .data$nFeature_RNA,
        colour = .data$classification
      )
    ) +
      ggplot2::geom_point(
        alpha = 0.34,
        size = 0.55
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(
        values = c(
          "Singlet" = "#4C78A8",
          "Doublet" = "#E45756"
        ),
        name = NULL
      ) +
      ggplot2::labs(
        title = "Library-complexity context",
        x = "RNA counts per cell (log10 scale)",
        y = "Detected genes per cell (log10 scale)"
      ) +
      ggplot2::theme_minimal(
        base_size = 10.5
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    p_doublet_diagnostic <-
      (p_doublet_score | p_doublet_complexity) +
      patchwork::plot_annotation(
        title = "Transcriptome-based doublet diagnostics",
        subtitle = paste0(
          "Complexity scatter is sampled for readability; classification uses all cells across ",
          dplyr::n_distinct(
            doublet_cell_diagnostics$sample_id
          ),
          " samples"
        ),
        caption = "High RNA content alone is not sufficient evidence of a doublet."
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "doublet_score_and_complexity.png"
      ),
      p_doublet_diagnostic,
      width = 12.5,
      height = 6.5,
      dpi = 300
    )


    if (REMOVE_DOUBLETS) {

      keep <- colnames(
        seurat_qc
      )[
        seurat_qc$scDblFinder.class ==
          "singlet"
      ]

      seurat_singlet <- subset(
        seurat_qc,
        cells = keep
      )

    } else {

      seurat_singlet <- seurat_qc
    }
  }


  save_checkpoint(
    seurat_singlet,
    "M08_singlets"
  )

  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M09. NORMALIZATION AND HVG SELECTION / 标准化与高变基因筛选
# ============================================================
#
# 【Module scope / 模块范围】
# This module separates library-size effects from relative expression, generates
# an exploratory normalized representation and selects genes with excess biological
# variation for dimensionality reduction. Raw RNA counts remain unchanged.
# 本模块校正文库大小差异，生成用于探索性分析的标准化表达表示，并筛选具有额外
# 生物学变异的基因用于降维；RNA原始计数保持不变。
#
# 【Quality gate / 质量门槛】
# Size-factor distributions, mean-variance trends, HVG counts and overlap with
# technical gene families should be reviewed. A normalization warning must be
# interpreted in the context of Seurat v5 layers and confirmed not to overwrite
# the raw count layer.
# 应检查size factor分布、均值—方差趋势、HVG数量及其与技术性基因家族的重叠。
# 对标准化warning应结合Seurat v5 layer结构判断，并确认未覆盖原始count layer。
#
# 【Interpretation boundary / 解读边界】
# HVGs define the features used to discover major cell-state structure; they are
# not a list of statistically significant condition-associated genes.
# HVG定义用于识别主要细胞状态结构的特征集合，并不代表与condition显著相关的基因。
#
# 【OSCA对应 / OSCA connection】
#
# OSCA basic - Normalization + Feature selection.
#
#
# 【默认路线 / Default route: OSCA_SCRAN】
#
# quickCluster()
#
# 先将表达相似的细胞进行粗分组。
# Performs coarse clustering before size-factor estimation.
#
#
# computeSumFactors()
#
# 通过pooling/deconvolution估计每个细胞的size factor，
# 对稀疏scRNA-seq数据通常比简单总count缩放更稳健。
#
# Estimates cell-specific size factors by pooling/deconvolution,
# addressing sparsity in single-cell count data.
#
#
# normalizeRnaCounts.se()
#
# 根据前面估计的size factor对counts进行标准化，
# 并生成用于后续探索性分析的log-normalized expression。
#
# Normalizes counts using the estimated size factors and
# generates log-normalized expression values.
#
#
# modelGeneVariances() / modelGeneVar() compatibility
#
# 对log-normalized expression的mean-variance关系进行建模，
# 将总方差分解为technical trend和biological component。
#
# Models the mean-variance relationship and estimates
# excess biological variation. Current Bioconductor releases use
# scrapper::modelGeneVariances(); the older scran interface remains a fallback.
#
#
# block = sample_id
#
# 多样本数据中，将sample/library作为blocking factor，
# 尽量避免某个样本整体表达差异主导HVG筛选。
#
# In multi-sample datasets, sample identity is used as a
# blocking factor to reduce the influence of sample-level
# shifts on HVG selection.
#
#
# chooseHighlyVariableGenes() / getTopHVGs() compatibility
#
# 根据biological component选择后续PCA和聚类最有信息的基因。
#
# Selects informative genes for downstream PCA and clustering. The current
# scrapper interface is used when available, with an older scran fallback.
#
#
# 【原始counts是否丢失？ / Are raw counts lost?】
#
# 不会。
#
# RNA raw counts继续保存在RNA assay的counts layer中。
# log-normalized expression用于PCA、聚类和可视化等探索性分析。
#
# No.
#
# Raw RNA counts remain available in the RNA assay.
# Log-normalized expression is used for exploratory analyses
# such as PCA, clustering and visualization.
#
#
# 【为什么需要兼容Seurat v5？ / Why Seurat v5 handling?】
#
# Seurat v5在merge多个样本后可能产生多个counts layers，例如：
#
# counts.Sample1
# counts.Sample2
# counts.Sample3
#
# OSCA/Bioconductor分析通常需要一个统一的count matrix。
#
# 因此在转换为SingleCellExperiment之前，
# 使用前面定义的prepare_for_sce()统一RNA layers。
#
#
# 【Output / 输出】
#
# OSCA_gene_variance.csv
#   OSCA/SCRAN路线下每个基因的variance modelling结果
#
# variable_genes.csv
#   最终选择的HVG
#
# HVG_mean_variance_trend.png
#   全部基因的均值—方差关系、technical trend及最终HVG
#
# normalization_size_factors_by_sample.png（OSCA_SCRAN路线）
#   各样本deconvolution size factor分布
#
# HVG_biological_component_top.png（OSCA_SCRAN路线）
#   biological variance component最高的代表性HVG
#
# OSCA_top_HVG_biological_variance.csv（OSCA_SCRAN路线）
#   带biological variance与技术基因类别的完整HVG审计表
#
# M09_normalized_HVG.rds
#   标准化并完成HVG选择的Seurat对象
#
#
# 【Check after running / 运行后检查】
#
# 1. 检查size factor是否存在明显异常。
#
# 2. 检查HVG modelling使用的sample block数量是否合理。
#
# 3. 检查top HVGs是否被mt-/Rpl/Rps等
#    technical/housekeeping signals完全支配。
#
# 4. 确认RNA raw counts仍然存在。
#
# ============================================================


if (RUN_M09_NORMALIZE_HVG) {

  MODULE <- "M09_normalize_HVG"

  MDIR <- module_dir(
    MODULE
  )


  # ==========================================================
  # 1. LOAD M08 OUTPUT
  # 1. 读取M08输出
  # ==========================================================

  seurat_singlet <- load_checkpoint(
    "M08_singlets"
  )


  DefaultAssay(
    seurat_singlet
  ) <- "RNA"


  # ==========================================================
  # 2. OSCA / SCRAN NORMALIZATION
  # 2. OSCA / SCRAN标准化路线
  # ==========================================================

  if (
    NORMALIZATION_METHOD ==
    "OSCA_SCRAN"
  ) {


    # --------------------------------------------------------
    # 2.1 Prepare Seurat v5 RNA assay
    # 2.1 整理Seurat v5 RNA assay
    # --------------------------------------------------------
    #
    # A merged Seurat v5 object may contain multiple
    # sample-specific counts layers.
    #
    # merge后的Seurat v5对象可能包含：
    #
    # counts.Sample1
    # counts.Sample2
    # ...
    #
    # prepare_for_sce() joins these layers before conversion
    # to SingleCellExperiment.
    #
    # prepare_for_sce()负责在转换SCE之前统一这些layers。
    # --------------------------------------------------------

    seurat_singlet <- prepare_for_sce(
      seurat_singlet,
      assay = "RNA"
    )


    # --------------------------------------------------------
    # 2.2 Convert Seurat to SingleCellExperiment
    # 2.2 转换为SingleCellExperiment
    # --------------------------------------------------------

    sce_norm <- as.SingleCellExperiment(
      seurat_singlet,
      assay = "RNA"
    )


    # --------------------------------------------------------
    # 2.3 Basic SCE validation
    # 2.3 SCE基本检查
    # --------------------------------------------------------

    if (
      !"counts" %in%
      SummarizedExperiment::assayNames(
        sce_norm
      )
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "RNA counts assay is missing after conversion to SCE. ",
          "转换为SCE后没有找到RNA counts assay。"
        )
      )
    }


    if (ncol(sce_norm) == 0L) {

      safe_stop_module(
        MODULE,
        "No cells remain before normalization / 标准化前已经没有cell。"
      )
    }


    cat(
      "\n============================================================\n",
      "OSCA/SCRAN NORMALIZATION\n",
      "OSCA/SCRAN 标准化\n",
      "============================================================\n",
      "Genes / 基因数: ",
      nrow(sce_norm),
      "\n",
      "Cells / cell数: ",
      ncol(sce_norm),
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # 2.4 Coarse clustering
    # 2.4 粗聚类
    # --------------------------------------------------------
    #
    # quickCluster() groups cells with broadly similar
    # expression profiles before deconvolution.
    #
    # quickCluster()首先将表达特征大致相似的cells
    # 进行粗分组。
    #
    # This improves pooling-based size-factor estimation
    # in heterogeneous scRNA-seq datasets.
    #
    # 对于异质性较强的单细胞数据，
    # 这样有助于后面的pooling/deconvolution。
    # --------------------------------------------------------

    set.seed(
      RANDOM_SEED
    )


    qc_cluster <- scran::quickCluster(
      sce_norm
    )


    if (
      length(qc_cluster) !=
      ncol(sce_norm)
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "quickCluster result length does not match cell number. ",
          "quickCluster结果长度与cell数量不一致。 ",
          "clusters = ",
          length(qc_cluster),
          "; cells = ",
          ncol(sce_norm)
        )
      )
    }


    cat(
      "\nquickCluster completed / quickCluster完成\n",
      "Coarse clusters / 粗分组数: ",
      length(unique(qc_cluster)),
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # 2.5 Estimate deconvolution size factors
    # 2.5 估计deconvolution size factors
    # --------------------------------------------------------
    #
    # computeSumFactors() estimates cell-specific size factors
    # by pooling cells and deconvolving the pooled factors.
    #
    # computeSumFactors()通过cell pooling和deconvolution
    # 估计每个cell的size factor。
    #
    # Raw counts are NOT overwritten.
    #
    # 原始counts不会被覆盖。
    # --------------------------------------------------------

    sce_norm <- scran::computeSumFactors(
      sce_norm,
      cluster = qc_cluster
    )


    sf <- SingleCellExperiment::sizeFactors(
      sce_norm
    )


    if (
      length(sf) !=
      ncol(sce_norm)
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "Size-factor length does not match cell number. ",
          "size factor数量与cell数量不一致。"
        )
      )
    }


    if (
      anyNA(sf) ||
      any(!is.finite(sf)) ||
      any(sf <= 0)
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "Invalid size factors detected. ",
          "检测到NA、非有限值或<=0的size factor。"
        )
      )
    }


    cat(
      "\nSize factors estimated / size factor估计完成\n"
    )


    print(
      summary(sf)
    )


    # --------------------------------------------------------
    # 2.6 Log-normalize RNA counts
    # 2.6 对RNA counts进行log normalization
    # --------------------------------------------------------
    #
    # The current scrapper API provides
    # normalizeRnaCounts.se() for
    # SummarizedExperiment/SingleCellExperiment input.
    #
    # 当前scrapper版本中，
    # normalizeRnaCounts.se()用于处理
    # SummarizedExperiment/SingleCellExperiment对象。
    #
    # The size factors estimated above are already stored
    # inside sce_norm.
    #
    # 前面估计的size factor已经储存在sce_norm中。
    #
    # The normalized expression should be written to
    # a "logcounts" assay.
    #
    # 标准化后的表达矩阵应保存为logcounts assay。
    # --------------------------------------------------------

    sce_norm <- scrapper::normalizeRnaCounts.se(
      sce_norm
    )


    if (
      !"logcounts" %in%
      SummarizedExperiment::assayNames(
        sce_norm
      )
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "Normalization did not generate a logcounts assay. ",
          "标准化后没有生成logcounts assay。"
        )
      )
    }


    cat(
      "\nNormalization completed / 标准化完成\n",
      "SCE assays / 当前SCE assays: ",
      paste(
        SummarizedExperiment::assayNames(
          sce_norm
        ),
        collapse = ", "
      ),
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # 2.7 Prepare sample-aware blocking factor
    # 2.7 准备sample-aware blocking factor
    # --------------------------------------------------------
    #
    # Gene-variance modelling can operate separately
    # across blocks and combine the results.
    #
    # Gene-variance modelling可以利用block减少样本整体差异
    # 对HVG筛选的影响。
    #
    # IMPORTANT:
    #
    # The block vector must contain exactly ONE value
    # for every column/cell in sce_norm.
    #
    # block必须与sce_norm中的cell严格一一对应。
    #
    # Therefore sample identity is extracted DIRECTLY
    # from colData(sce_norm), rather than separately
    # from the Seurat object.
    #
    # 因此这里直接从真正进入gene-variance model的
    # sce_norm colData中读取sample_id。
    #
    # This avoids the previous problem where Seurat [[ ]]
    # could return a one-column data.frame rather than
    # a simple vector.
    #
    # 这样可以避免Seurat [[ ]]返回单列data.frame，
    # 导致block长度判断错误的问题。
    # --------------------------------------------------------

    sce_coldata <- SummarizedExperiment::colData(
      sce_norm
    )


    if (
      !SAMPLE_ID_COL %in%
      colnames(sce_coldata)
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "Sample ID column not found in SCE metadata / ",
          "SCE metadata中没有找到sample ID列: ",
          SAMPLE_ID_COL
        )
      )
    }


    block <- as.character(
      sce_coldata[[SAMPLE_ID_COL]]
    )


    # --------------------------------------------------------
    # 2.8 Strict block validation
    # 2.8 严格检查block
    # --------------------------------------------------------

    if (
      length(block) !=
      ncol(sce_norm)
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "HVG blocking vector length does not match ",
          "the number of cells in sce_norm.\n",
          "HVG block长度与sce_norm cell数量不一致。\n",
          "block = ",
          length(block),
          "; cells = ",
          ncol(sce_norm)
        )
      )
    }


    if (anyNA(block)) {

      safe_stop_module(
        MODULE,
        paste0(
          "Sample ID contains NA values during HVG modelling / ",
          "HVG建模时sample ID存在NA: ",
          SAMPLE_ID_COL
        )
      )
    }


    if (any(trimws(block) == "")) {

      safe_stop_module(
        MODULE,
        paste0(
          "Sample ID contains empty values during HVG modelling / ",
          "HVG建模时sample ID存在空值: ",
          SAMPLE_ID_COL
        )
      )
    }


    block <- factor(
      block
    )


    cat(
      "\n============================================================\n",
      "HVG VARIANCE MODELLING BLOCKS\n",
      "HVG 方差建模分组\n",
      "============================================================\n",
      "Cells / cell数: ",
      ncol(sce_norm),
      "\n",
      "Blocks / sample数: ",
      nlevels(block),
      "\n\n",
      sep = ""
    )


    print(
      table(block)
    )


    # --------------------------------------------------------
    # 2.9 Model gene variance
    # 2.9 建模基因方差
    # --------------------------------------------------------
    #
    # THIS STEP CREATES gene_var.
    #
    # 这一步非常重要：
    # gene_var就是在这里正式生成的。
    #
    # The previous version accidentally removed this section,
    # while later code still attempted to use gene_var.
    #
    # 前一个版本中这一段被意外删除，
    # 但后面仍然使用gene_var，
    # 因此出现：
    #
    # object 'gene_var' not found
    #
    # Current Bioconductor releases provide modelGeneVariances() and
    # chooseHighlyVariableGenes() in scrapper. Older releases used
    # scran::modelGeneVar() and scran::getTopHVGs(). The modern route is
    # preferred when available; a legacy fallback keeps the template usable
    # in older, otherwise compatible Bioconductor installations.
    #
    # 当前Bioconductor推荐使用scrapper::modelGeneVariances()和
    # chooseHighlyVariableGenes()。较旧版本使用scran::modelGeneVar()和
    # getTopHVGs()。代码优先采用新接口，并保留旧版本fallback，从而兼容
    # 不同Bioconductor版本而不改变“sample-aware variance modelling”原则。
    # --------------------------------------------------------

    scrapper_exports <- getNamespaceExports(
      "scrapper"
    )

    use_modern_HVG_API <- all(
      c(
        "modelGeneVariances",
        "chooseHighlyVariableGenes"
      ) %in% scrapper_exports
    )

    if (use_modern_HVG_API) {

      modern_gene_var <- scrapper::modelGeneVariances(
        SummarizedExperiment::assay(
          sce_norm,
          "logcounts"
        ),
        block = block
      )

      modern_gene_statistics <- as.data.frame(
        modern_gene_var$statistics
      )

      # Harmonize current column names with the long-standing OSCA terminology
      # used elsewhere in this template and in exported audit tables.
      # 将新接口列名统一为母版长期使用的OSCA术语：
      # means -> mean; variances -> total; fitted -> tech;
      # residuals -> bio。
      gene_var <- data.frame(
        mean = modern_gene_statistics$means,
        total = modern_gene_statistics$variances,
        tech = modern_gene_statistics$fitted,
        bio = modern_gene_statistics$residuals,
        row.names = rownames(
          modern_gene_statistics
        ),
        stringsAsFactors = FALSE
      )

    } else {

      gene_var <- scran::modelGeneVar(
        sce_norm,
        block = block
      )
    }


    if (
      is.null(gene_var) ||
      nrow(gene_var) == 0L
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "modelGeneVar returned no genes / ",
          "modelGeneVar没有返回有效基因。"
        )
      )
    }


    cat(
      "\nGene variance modelling completed / ",
      "基因方差建模完成\n",
      "Genes modelled / 完成建模的基因数: ",
      nrow(gene_var),
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # 2.10 Determine the actual number of HVGs
    # 2.10 确定实际选择的HVG数量
    # --------------------------------------------------------
    #
    # N_HVG is the requested maximum.
    #
    # N_HVG是用户设置的目标HVG数量。
    #
    # If fewer genes are available, all available genes
    # up to that number are used.
    #
    # 如果可用基因少于N_HVG，则自动使用实际可用数量。
    # --------------------------------------------------------

    n_hvg_actual <- min(
      as.integer(N_HVG),
      nrow(gene_var)
    )


    if (n_hvg_actual < 2L) {

      safe_stop_module(
        MODULE,
        paste0(
          "Fewer than two genes are available for HVG selection. ",
          "可用于HVG筛选的基因少于2个。"
        )
      )
    }


    # --------------------------------------------------------
    # 2.11 Select top HVGs
    # 2.11 选择top HVGs
    # --------------------------------------------------------

    if (use_modern_HVG_API) {

      hvg_index <- scrapper::chooseHighlyVariableGenes(
        stats = as.numeric(
          gene_var$bio
        ),
        top = n_hvg_actual,
        larger = TRUE,
        keep.ties = FALSE,
        bound = 0
      )

      hvg <- rownames(
        gene_var
      )[
        hvg_index
      ]

    } else {

      hvg <- scran::getTopHVGs(
        gene_var,
        n = n_hvg_actual
      )
    }


    hvg <- unique(
      as.character(hvg)
    )


    # Keep only genes that are actually present
    # in the Seurat RNA assay.
    #
    # 再次确认HVG确实存在于Seurat RNA assay中。

    hvg <- intersect(
      hvg,
      rownames(
        seurat_singlet[["RNA"]]
      )
    )


    if (length(hvg) < 2L) {

      safe_stop_module(
        MODULE,
        paste0(
          "Too few HVGs remain after matching to the Seurat object. ",
          "与Seurat对象匹配后剩余HVG过少。"
        )
      )
    }


    cat(
      "\nHVG selection completed / HVG筛选完成\n",
      "HVGs retained / 保留HVG数: ",
      length(hvg),
      "\n",
      sep = ""
    )


    # --------------------------------------------------------
    # 2.12 Extract log-normalized expression matrix
    # 2.12 提取log-normalized expression矩阵
    # --------------------------------------------------------

    logcounts_matrix <- SummarizedExperiment::assay(
      sce_norm,
      "logcounts"
    )

    # Materialize the delayed matrix before storing it in Seurat.
    # 在写入Seurat之前，将延迟矩阵显式转换为可序列化稀疏矩阵。
    logcounts_matrix <- methods::as(
      logcounts_matrix,
      "dgCMatrix"
    )


    # --------------------------------------------------------
    # 2.13 Confirm exact cell correspondence
    # 2.13 确认SCE和Seurat cell顺序完全一致
    # --------------------------------------------------------
    #
    # Before writing logcounts back to Seurat,
    # the cell names and their order must match.
    #
    # 在把logcounts写回Seurat之前，
    # 必须确保两个对象中的cell名称和顺序完全一致。
    # --------------------------------------------------------

    if (
      !identical(
        colnames(logcounts_matrix),
        colnames(seurat_singlet)
      )
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "Cell names/order differ between SCE logcounts and Seurat. ",
          "SCE logcounts与Seurat的cell名称或顺序不一致。"
        )
      )
    }


    if (
      !identical(
        rownames(logcounts_matrix),
        rownames(seurat_singlet[["RNA"]])
      )
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "Gene names/order differ between SCE logcounts and Seurat. ",
          "SCE logcounts与Seurat RNA assay的基因名称或顺序不一致。"
        )
      )
    }


    # --------------------------------------------------------
    # 2.14 Write normalized data back to Seurat
    # 2.14 将标准化表达量写回Seurat
    # --------------------------------------------------------
    #
    # Seurat v5 uses layers.
    #
    # Seurat v5使用layer储存不同形式的数据：
    #
    # counts
    #   raw counts
    #
    # data
    #   normalized expression
    #
    # scale.data
    #   scaled expression used for PCA
    #
    # We first try the Seurat v5 layer interface.
    #
    # 这里优先使用Seurat v5的layer接口。
    #
    # The slot fallback is retained only for compatibility
    # with older SeuratObject versions.
    #
    # slot仅作为旧版本兼容fallback。
    # --------------------------------------------------------

    set_data_ok <- FALSE


    try({

      seurat_singlet <- SeuratObject::SetAssayData(
        seurat_singlet,
        assay = "RNA",
        layer = "data",
        new.data = logcounts_matrix
      )

      set_data_ok <- TRUE

    }, silent = TRUE)


    if (!set_data_ok) {

      seurat_singlet <- SeuratObject::SetAssayData(
        seurat_singlet,
        assay = "RNA",
        slot = "data",
        new.data = logcounts_matrix
      )
    }


    # --------------------------------------------------------
    # 2.15 Register HVGs in Seurat
    # 2.15 将HVG写入Seurat
    # --------------------------------------------------------

    VariableFeatures(
      seurat_singlet
    ) <- hvg


    # --------------------------------------------------------
    # 2.16 Decide whether to regress percent.mt
    # 2.16 判断是否回归percent.mt
    # --------------------------------------------------------
    #
    # Regression is optional.
    #
    # 是否回归percent.mt由全局参数REGRESS_PERCENT_MT决定。
    #
    # It should not automatically be performed simply because
    # percent.mt exists.
    #
    # 不能因为存在percent.mt就自动进行回归。
    # --------------------------------------------------------

    vars_to_regress <- if (
      REGRESS_PERCENT_MT
    ) {

      if (
        !"percent.mt" %in%
        colnames(seurat_singlet[[]])
      ) {

        safe_stop_module(
          MODULE,
          paste0(
            "REGRESS_PERCENT_MT=TRUE but percent.mt is missing. ",
            "要求回归percent.mt，但metadata中没有percent.mt。"
          )
        )
      }

      "percent.mt"

    } else {

      NULL
    }


    # --------------------------------------------------------
    # 2.17 Scale selected HVGs
    # 2.17 对HVG进行ScaleData
    # --------------------------------------------------------
    #
    # Only selected HVGs are scaled here because they are
    # the features used for downstream PCA.
    #
    # 这里只scale选中的HVG，
    # 因为后面的PCA主要使用这些基因。
    # --------------------------------------------------------

    seurat_singlet <- ScaleData(
      seurat_singlet,
      features = hvg,
      vars.to.regress = vars_to_regress,
      verbose = TRUE
    )


    ANALYSIS_ASSAY <- "RNA"


    # --------------------------------------------------------
    # 2.18 Save OSCA gene-variance results
    # 2.18 保存OSCA gene variance结果
    # --------------------------------------------------------

    gene_var_df <- as.data.frame(
      gene_var
    )


    gene_var_df$gene <- rownames(
      gene_var_df
    )


    # Put gene name first for easier inspection.
    # 将gene放在第一列，方便人工查看。

    gene_var_df <- gene_var_df[
      ,
      c(
        "gene",
        setdiff(
          colnames(gene_var_df),
          "gene"
        )
      ),
      drop = FALSE
    ]


    data.table::fwrite(
      gene_var_df,
      file.path(
        MDIR,
        "tables",
        "OSCA_gene_variance.csv"
      )
    )


    # ==========================================================
    # 3. SEURAT LOGNORMALIZE ROUTE
    # 3. Seurat LogNormalize路线
    # ==========================================================

  } else if (
    NORMALIZATION_METHOD ==
    "LOGNORMALIZE"
  ) {


    # --------------------------------------------------------
    # Standard Seurat normalization
    # 标准Seurat LogNormalize
    # --------------------------------------------------------

    seurat_singlet <- NormalizeData(
      seurat_singlet,
      normalization.method = "LogNormalize",
      scale.factor = 10000
    )


    # --------------------------------------------------------
    # Select variable features using Seurat VST
    # 使用Seurat VST筛选HVG
    # --------------------------------------------------------

    seurat_singlet <- FindVariableFeatures(
      seurat_singlet,
      selection.method = "vst",
      nfeatures = min(
        as.integer(N_HVG),
        nrow(seurat_singlet)
      )
    )


    # --------------------------------------------------------
    # Optional mitochondrial regression
    # 可选percent.mt回归
    # --------------------------------------------------------

    vars_to_regress <- if (
      REGRESS_PERCENT_MT
    ) {

      if (
        !"percent.mt" %in%
        colnames(seurat_singlet[[]])
      ) {

        safe_stop_module(
          MODULE,
          paste0(
            "REGRESS_PERCENT_MT=TRUE but percent.mt is missing. ",
            "要求回归percent.mt，但metadata中没有percent.mt。"
          )
        )
      }

      "percent.mt"

    } else {

      NULL
    }


    seurat_singlet <- ScaleData(
      seurat_singlet,
      features = VariableFeatures(
        seurat_singlet
      ),
      vars.to.regress = vars_to_regress
    )


    ANALYSIS_ASSAY <- "RNA"


    # ==========================================================
    # 4. SCT ROUTE
    # 4. SCTransform路线
    # ==========================================================

  } else if (
    NORMALIZATION_METHOD ==
    "SCT"
  ) {


    vars_to_regress <- if (
      REGRESS_PERCENT_MT
    ) {

      if (
        !"percent.mt" %in%
        colnames(seurat_singlet[[]])
      ) {

        safe_stop_module(
          MODULE,
          paste0(
            "REGRESS_PERCENT_MT=TRUE but percent.mt is missing. ",
            "要求回归percent.mt，但metadata中没有percent.mt。"
          )
        )
      }

      "percent.mt"

    } else {

      NULL
    }


    seurat_singlet <- SCTransform(
      seurat_singlet,
      assay = "RNA",
      new.assay.name = "SCT",
      vars.to.regress = vars_to_regress,
      variable.features.n = min(
        as.integer(N_HVG),
        nrow(seurat_singlet)
      ),
      vst.flavor = "v2",
      seed.use = RANDOM_SEED
    )


    DefaultAssay(
      seurat_singlet
    ) <- "SCT"


    ANALYSIS_ASSAY <- "SCT"


    # ==========================================================
    # 5. UNSUPPORTED METHOD
    # 5. 不支持的标准化方法
    # ==========================================================

  } else {

    safe_stop_module(
      MODULE,
      paste0(
        "Unsupported NORMALIZATION_METHOD / ",
        "不支持的NORMALIZATION_METHOD: ",
        NORMALIZATION_METHOD
      )
    )
  }


  # ==========================================================
  # 6. RECORD ANALYSIS ASSAY
  # 6. 记录后续分析使用的assay
  # ==========================================================
  #
  # M10 and later modules can retrieve this attribute
  # instead of guessing which assay should be used.
  #
  # M10及后续模块可以直接读取该attribute，
  # 而不需要重新猜测应该使用RNA还是SCT。
  # ==========================================================

  attr(
    seurat_singlet,
    "ANALYSIS_ASSAY"
  ) <- ANALYSIS_ASSAY


  # ==========================================================
  # 7. SAVE VARIABLE GENES
  # 7. 保存HVG列表
  # ==========================================================

  final_hvg <- VariableFeatures(
    seurat_singlet
  )


  if (length(final_hvg) == 0L) {

    safe_stop_module(
      MODULE,
      paste0(
        "No variable genes were retained / ",
        "没有获得任何variable genes。"
      )
    )
  }


  variable_genes <- data.frame(
    rank = seq_along(
      final_hvg
    ),
    gene = final_hvg,
    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    variable_genes,
    file.path(
      MDIR,
      "tables",
      "variable_genes.csv"
    )
  )


  # ==========================================================
  # 8. NORMALIZATION AND HVG DIAGNOSTIC FIGURES
  # 8. 标准化与HVG诊断图片
  # ==========================================================
  #
  # 图中可见文字使用英文，避免不同PNG设备缺少中文字形时出现乱码；
  # 方法说明和代码注释继续保持详细中英双语。所有图片仅汇报本模块
  # 已计算的统计量，不改变normalization、HVG筛选或后续PCA输入。
  #
  # Visible plot text is English-only for portable rendering. These figures
  # report quantities already calculated by M09 and do not alter normalization,
  # feature selection or the genes passed to PCA.
  # ==========================================================

  if (
    NORMALIZATION_METHOD ==
      "OSCA_SCRAN"
  ) {

    # --------------------------------------------------------
    # 8.1 MEAN-VARIANCE TREND
    # 8.1 均值—方差趋势
    # --------------------------------------------------------
    #
    # total表示每个基因观测到的总方差，tech表示modelGeneVar拟合的
    # technical trend，bio = total - tech表示超出技术趋势的方差成分。
    # 蓝色点是最终进入PCA特征空间的HVG；红线不是线性回归，而是
    # scran模型估计的technical trend。
    #
    # total is the observed variance, tech is the fitted technical trend and
    # bio is the excess component. The red curve is the scran model output,
    # not an independently fitted straight regression line.
    # --------------------------------------------------------

    required_gene_var_columns <- c(
      "gene",
      "mean",
      "total",
      "tech",
      "bio"
    )

    missing_gene_var_columns <- setdiff(
      required_gene_var_columns,
      colnames(gene_var_df)
    )

    if (
      length(
        missing_gene_var_columns
      ) > 0L
    ) {
      safe_stop_module(
        MODULE,
        paste0(
          "OSCA variance table lacks columns required for plotting: ",
          paste(
            missing_gene_var_columns,
            collapse = ", "
          ),
          " / OSCA variance结果缺少绘图所需列。"
        )
      )
    }

    gene_variance_plot_data <- gene_var_df |>
      dplyr::mutate(
        selected_HVG = ifelse(
          .data$gene %in%
            final_hvg,
          "Selected HVG",
          "Other gene"
        )
      ) |>
      dplyr::filter(
        is.finite(
          .data$mean
        ),
        is.finite(
          .data$total
        ),
        is.finite(
          .data$tech
        ),
        .data$mean > 0,
        .data$total > 0,
        .data$tech > 0
      )

    labelled_HVGs <- gene_variance_plot_data |>
      dplyr::filter(
        .data$selected_HVG ==
          "Selected HVG",
        is.finite(
          .data$bio
        )
      ) |>
      dplyr::arrange(
        dplyr::desc(
          .data$bio
        )
      ) |>
      dplyr::slice_head(
        n = min(
          as.integer(
            HVG_DIAGNOSTIC_LABEL_TOP
          ),
          nrow(
            gene_variance_plot_data
          )
        )
      )

    technical_trend_data <- gene_variance_plot_data |>
      dplyr::arrange(
        .data$mean
      )

    p_hvg_mean_variance <- ggplot2::ggplot(
      gene_variance_plot_data,
      ggplot2::aes(
        x = .data$mean,
        y = .data$total
      )
    ) +
      ggplot2::geom_point(
        ggplot2::aes(
          colour = .data$selected_HVG
        ),
        alpha = 0.44,
        size = 0.62
      ) +
      ggplot2::geom_line(
        data = technical_trend_data,
        ggplot2::aes(
          x = .data$mean,
          y = .data$tech
        ),
        inherit.aes = FALSE,
        colour = "#D94841",
        linewidth = 0.88,
        alpha = 0.95
      ) +
      ggrepel::geom_text_repel(
        data = labelled_HVGs,
        ggplot2::aes(
          label = .data$gene
        ),
        seed = RANDOM_SEED,
        size = 3,
        colour = "#17365D",
        box.padding = 0.35,
        point.padding = 0.18,
        min.segment.length = 0,
        max.overlaps = Inf,
        segment.colour = "#8C8C8C",
        show.legend = FALSE
      ) +
      ggplot2::scale_x_log10() +
      ggplot2::scale_y_log10() +
      ggplot2::scale_colour_manual(
        values = c(
          "Other gene" = "#B8B8B8",
          "Selected HVG" = "#2F6BCE"
        ),
        breaks = c(
          "Selected HVG",
          "Other gene"
        ),
        name = NULL
      ) +
      ggplot2::labs(
        title = "HVG mean-variance trend",
        subtitle = paste0(
          "Observed gene variance versus mean expression; red = fitted technical trend; selected HVGs = ",
          length(
            final_hvg
          )
        ),
        x = "Mean log-expression (log10 scale)",
        y = "Observed variance (log10 scale)",
        caption = "HVGs reflect excess biological variability and are not condition-specific differential-expression results."
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "HVG_mean_variance_trend.png"
      ),
      p_hvg_mean_variance,
      width = 9.5,
      height = 7.2,
      dpi = 300
    )


    # --------------------------------------------------------
    # 8.2 SIZE-FACTOR DISTRIBUTIONS BY SAMPLE
    # 8.2 各样本size factor分布
    # --------------------------------------------------------
    #
    # 使用violin显示完整分布，内部boxplot显示中位数和四分位数。
    # y轴采用log10尺度，便于识别跨数量级异常。样本间位置差异可能来自
    # 测序深度或细胞组成，不应自动解释为batch effect。
    # --------------------------------------------------------

    size_factor_plot_data <- data.frame(
      cell = colnames(
        sce_norm
      ),
      sample_id = as.character(
        block
      ),
      condition = if (
        CONDITION_COL %in%
          colnames(sce_coldata)
      ) {
        as.character(
          sce_coldata[[CONDITION_COL]]
        )
      } else {
        "Unspecified"
      },
      size_factor = as.numeric(
        sf
      ),
      stringsAsFactors = FALSE
    ) |>
      dplyr::filter(
        is.finite(
          .data$size_factor
        ),
        .data$size_factor > 0
      )

    p_size_factors <- ggplot2::ggplot(
      size_factor_plot_data,
      ggplot2::aes(
        x = .data$sample_id,
        y = .data$size_factor,
        fill = .data$condition
      )
    ) +
      ggplot2::geom_violin(
        scale = "width",
        trim = TRUE,
        alpha = 0.72,
        colour = NA
      ) +
      ggplot2::geom_boxplot(
        width = 0.14,
        outlier.shape = NA,
        fill = "white",
        colour = "#333333",
        linewidth = 0.34
      ) +
      ggplot2::scale_y_log10() +
      ggplot2::labs(
        title = "Deconvolution size factors by sample",
        subtitle = "Violin distributions with median and interquartile range; OSCA/scran normalization",
        x = NULL,
        y = "Cell-specific size factor (log10 scale)",
        fill = "Condition",
        caption = "Extreme values should be reviewed with RNA counts, detected genes and sample-level QC."
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "normalization_size_factors_by_sample.png"
      ),
      p_size_factors,
      width = max(
        9,
        0.65 *
          dplyr::n_distinct(
            size_factor_plot_data$sample_id
          ) +
          4
      ),
      height = 6.5,
      dpi = 300
    )


    # --------------------------------------------------------
    # 8.3 TOP BIOLOGICAL-VARIANCE HVGs
    # 8.3 biological variance最高的代表性HVG
    # --------------------------------------------------------
    #
    # 该图展示HVG中bio component最高的30个基因，并标记mitochondrial与
    # ribosomal类别，帮助判断HVG是否被技术/housekeeping信号支配。
    # 分类只用于审计，不会自动删除任何基因。
    # --------------------------------------------------------

    OSCA_top_HVG_table <- gene_var_df |>
      dplyr::filter(
        .data$gene %in%
          final_hvg
      ) |>
      dplyr::mutate(
        gene_category = dplyr::case_when(
          grepl(
            species_cfg$mito_pattern,
            .data$gene
          ) ~ "Mitochondrial",
          grepl(
            species_cfg$ribo_pattern,
            .data$gene
          ) ~ "Ribosomal",
          TRUE ~ "Other"
        )
      ) |>
      dplyr::arrange(
        dplyr::desc(
          .data$bio
        )
      ) |>
      dplyr::mutate(
        biological_variance_rank = dplyr::row_number()
      )

    data.table::fwrite(
      OSCA_top_HVG_table,
      file.path(
        MDIR,
        "tables",
        "OSCA_top_HVG_biological_variance.csv"
      )
    )

    top_biological_HVG_plot <- OSCA_top_HVG_table |>
      dplyr::filter(
        is.finite(
          .data$bio
        ),
        .data$bio > 0
      ) |>
      dplyr::slice_head(
        n = min(
          30L,
          nrow(
            OSCA_top_HVG_table
          )
        )
      ) |>
      dplyr::mutate(
        gene = factor(
          .data$gene,
          levels = rev(
            .data$gene
          )
        )
      )

    if (
      nrow(
        top_biological_HVG_plot
      ) > 0L
    ) {

      p_top_biological_HVG <- ggplot2::ggplot(
        top_biological_HVG_plot,
        ggplot2::aes(
          x = .data$bio,
          y = .data$gene,
          colour = .data$gene_category
        )
      ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0,
            xend = .data$bio,
            yend = .data$gene
          ),
          linewidth = 0.55,
          alpha = 0.72
        ) +
        ggplot2::geom_point(
          size = 2.8
        ) +
        ggplot2::scale_colour_manual(
          values = c(
            "Other" = "#2F6BCE",
            "Ribosomal" = "#F2A541",
            "Mitochondrial" = "#D94841"
          ),
          drop = FALSE,
          name = "Gene category"
        ) +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(
            mult = c(
              0,
              0.08
            )
          )
        ) +
        ggplot2::labs(
          title = "Leading HVGs by biological variance component",
          subtitle = "Top 30 selected HVGs ranked by variance above the fitted technical trend",
          x = "Estimated biological variance component",
          y = NULL,
          caption = "Mitochondrial and ribosomal labels are audit flags; genes are not automatically excluded."
        ) +
        ggplot2::theme_minimal(
          base_size = 10.5
        ) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(
            face = "bold"
          ),
          legend.position = "bottom"
        )

      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "HVG_biological_component_top.png"
        ),
        p_top_biological_HVG,
        width = 9,
        height = 9.5,
        dpi = 300
      )
    }

  } else {

    # LOGNORMALIZE和SCT路线没有OSCA的tech/bio分解，因此使用Seurat本身的
    # standardized variance图；不会把另一种模型的technical trend强行套入。
    # LOGNORMALIZE and SCT use Seurat's own feature-selection diagnostics because
    # an OSCA technical/biological variance decomposition is not available.
    p_hvg_mean_variance <- Seurat::VariableFeaturePlot(
      seurat_singlet,
      assay = ANALYSIS_ASSAY
    )

    label_features <- utils::head(
      final_hvg,
      n = min(
        as.integer(
          HVG_DIAGNOSTIC_LABEL_TOP
        ),
        length(
          final_hvg
        )
      )
    )

    if (length(label_features) > 0L) {
      p_hvg_mean_variance <- Seurat::LabelPoints(
        plot = p_hvg_mean_variance,
        points = label_features,
        repel = TRUE
      )
    }

    p_hvg_mean_variance <- p_hvg_mean_variance +
      ggplot2::labs(
        title = "HVG mean-variance diagnostic",
        subtitle = paste0(
          "Feature-selection method: ",
          NORMALIZATION_METHOD,
          "; selected HVGs = ",
          length(
            final_hvg
          )
        ),
        caption = "The plotted statistic follows the selected Seurat normalization route and is not differential expression."
      ) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "HVG_mean_variance_trend.png"
      ),
      p_hvg_mean_variance,
      width = 9.5,
      height = 7.2,
      dpi = 300
    )
  }


  # ==========================================================
  # 9. FINAL M09 SUMMARY
  # 9. M09最终汇总
  # ==========================================================

  cat(
    "\n============================================================\n",
    "M09 NORMALIZATION AND HVG SELECTION COMPLETE\n",
    "M09 标准化与HVG筛选完成\n",
    "============================================================\n",
    "Normalization method / 标准化方法: ",
    NORMALIZATION_METHOD,
    "\n",
    "Analysis assay / 后续分析assay: ",
    ANALYSIS_ASSAY,
    "\n",
    "Cells / cell数: ",
    ncol(seurat_singlet),
    "\n",
    "Genes / 基因数: ",
    nrow(seurat_singlet),
    "\n",
    "Variable genes / HVG数: ",
    length(final_hvg),
    "\n",
    "Raw RNA counts retained / 原始RNA counts保留: YES\n",
    "============================================================\n",
    sep = ""
  )


  # ==========================================================
  # 10. SAVE CHECKPOINT
  # 10. 保存checkpoint
  # ==========================================================

  save_checkpoint(
    seurat_singlet,
    "M09_normalized_HVG"
  )


  write_module_status(
    MODULE,
    "OK"
  )
}


# ============================================================
# M10. PCA, NEIGHBOR GRAPH, CLUSTERING AND UMAP / PCA、邻居图、聚类与UMAP
# ============================================================
#
# 【Module scope / 模块范围】
# PCA summarizes variation across the selected HVGs, the neighbor graph encodes
# local transcriptional similarity, graph clustering proposes candidate cell
# populations, and UMAP provides a two-dimensional visualization of that graph.
# PCA概括所选HVG中的主要变化，邻居图表示局部转录相似性，图聚类提出候选细胞群，
# UMAP则将该结构可视化为二维图形。
#
# 【Quality gate / 质量门槛】
# The number of PCs should be supported by elbow/variance diagnostics. Cluster
# resolution must be assessed with marker coherence, cluster size and resolution
# sensitivity, not selected solely because one UMAP appears visually attractive.
# PC数量应由elbow/方差诊断支持。聚类resolution需结合marker一致性、cluster大小和
# 分辨率敏感性评估，不能仅依据某一张UMAP的视觉效果选择。
#
# 【Interpretation boundary / 解读边界】
# Cluster numbers are arbitrary identifiers, UMAP axes have no direct biological
# units, and visual separation does not establish disease association or statistical
# significance. Formal case-control inference is deferred to replicate-level steps.
# cluster编号只是任意标识，UMAP坐标轴没有直接生物学单位，视觉分离也不能证明疾病
# 关联或统计显著性。正式病例—对照推断留待生物学重复层面的下游模块完成。
#
# 【OSCA对应 / OSCA connection】
#
# PCA用于把高维基因表达压缩到少量主要变化轴。
# PCA compresses high-dimensional expression into major axes
# of biological variation.
#
# 本模块使用M09已经确定的HVG进行：
#
# HVGs
#   ↓
# PCA
#   ↓
# KNN/SNN graph
#   ↓
# graph clustering
#   ↓
# UMAP
#
#
# 【重要原则 / Important principle】
#
# PCA / clustering / UMAP属于探索性分析。
#
# PCA, clustering and UMAP are exploratory analyses.
#
# UMAP上的距离、cluster大小或者Control/STZ之间的视觉差异
# 不能直接解释为统计显著性。
#
# Visual separation on UMAP is NOT a statistical test of
# condition-level differences.
#
#
# 【Input / 输入】
#
# M09_normalized_HVG.rds
#
#
# 【Output / 输出】
#
# PCA_elbow.png
# UMAP_overview.png
# UMAP_split_by_sample.png
# UMAP_cell_cycle.png（如果Phase metadata存在 / if Phase metadata is available）
# cluster_by_sample.csv
# M10_clustered.rds
#
#
# 【Check after running / 运行后检查】
#
# 1. cluster是否被单一样本垄断？
# 2. UMAP是否明显由sample/batch主导？
# 3. Control/Case是否与sample或batch完全混淆？
# 4. ElbowPlot是否支持当前N_PCS_USE？
#
# ============================================================


if (RUN_M10_PCA_CLUSTER_UMAP) {

  MODULE <- "M10_PCA_cluster_UMAP"
  MDIR <- module_dir(MODULE)


  # ==========================================================
  # 1. LOAD M09 OBJECT
  # 1. 读取M09对象
  # ==========================================================

  seurat_norm <- load_checkpoint(
    "M09_normalized_HVG"
  )


  # M09已经把实际分析assay记录在object attribute中。
  # M09 records the assay used for downstream analysis.

  ANALYSIS_ASSAY <- attr(
    seurat_norm,
    "ANALYSIS_ASSAY"
  )


  if (
    is.null(ANALYSIS_ASSAY) ||
    !ANALYSIS_ASSAY %in% SeuratObject::Assays(seurat_norm)
  ) {

    ANALYSIS_ASSAY <- DefaultAssay(
      seurat_norm
    )
  }


  DefaultAssay(
    seurat_norm
  ) <- ANALYSIS_ASSAY

  # ==========================================================
  # CELL-CYCLE DIAGNOSTIC
  # 细胞周期诊断
  # ==========================================================
  #
  # Cell-cycle scores are calculated as a diagnostic only.
  # 这里只计算cell-cycle score用于检查，不默认进行regression。
  #
  # Cell-cycle regression may remove biologically meaningful
  # disease-associated proliferative signals.
  #
  # 如果后续发现PCA/cluster明显被cell cycle主导，
  # 再单独考虑是否进行敏感性分析。
  # ==========================================================

  s_genes_present <- unname(
    Seurat::CaseMatch(
      search = Seurat::cc.genes.updated.2019$s.genes,
      match = rownames(seurat_norm)
    )
  )

  g2m_genes_present <- unname(
    Seurat::CaseMatch(
      search = Seurat::cc.genes.updated.2019$g2m.genes,
      match = rownames(seurat_norm)
    )
  )

  if (
    length(s_genes_present) >= 5L &&
    length(g2m_genes_present) >= 5L
  ) {

    seurat_norm <- CellCycleScoring(
      seurat_norm,
      s.features = s_genes_present,
      g2m.features = g2m_genes_present,
      set.ident = FALSE
    )


    p_cell_cycle <- VlnPlot(
      seurat_norm,
      features = c(
        "S.Score",
        "G2M.Score"
      ),
      group.by = SAMPLE_ID_COL,
      pt.size = 0
    )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "cell_cycle_scores_by_sample.png"
      ),
      p_cell_cycle,
      width = 12,
      height = 6,
      dpi = 300
    )


    cell_cycle_summary <- as.data.frame(
      table(
        sample = seurat_norm[[]][[SAMPLE_ID_COL]],
        phase = seurat_norm$Phase
      )
    )


    data.table::fwrite(
      cell_cycle_summary,
      file.path(
        MDIR,
        "tables",
        "cell_cycle_by_sample.csv"
      )
    )

  } else {

    warning(
      paste0(
        "Too few canonical cell-cycle genes were detected; ",
        "cell-cycle diagnostic was skipped."
      )
    )
  }

  # ==========================================================
  # 2. CHECK HVGs AND CELL NUMBER
  # 2. 检查HVG和cell数量
  # ==========================================================

  hvg_use <- VariableFeatures(
    seurat_norm
  )


  n_hvg_available <- length(
    hvg_use
  )

  n_cells_available <- ncol(
    seurat_norm
  )


  if (n_hvg_available < 2L) {

    safe_stop_module(
      MODULE,
      paste0(
        "Too few variable genes for PCA / ",
        "用于PCA的HVG少于2个。"
      )
    )
  }


  if (n_cells_available < 3L) {

    safe_stop_module(
      MODULE,
      paste0(
        "Too few cells for PCA/clustering / ",
        "用于PCA/聚类的cell数量过少。"
      )
    )
  }


  # ==========================================================
  # 3. DETERMINE SAFE PCA DIMENSIONS
  # 3. 自动确定安全的PC数量
  # ==========================================================

  N_PCS_ACTUAL <- min(
    as.integer(N_PCS_COMPUTE),
    n_hvg_available - 1L,
    n_cells_available - 1L
  )


  if (N_PCS_ACTUAL < 2L) {

    safe_stop_module(
      MODULE,
      "Too few cells/HVGs to compute PCA safely."
    )
  }


  N_PCS_USE_ACTUAL <- min(
    as.integer(N_PCS_USE),
    N_PCS_ACTUAL
  )


  dims_use <- seq_len(
    N_PCS_USE_ACTUAL
  )


  # Neighbor parameters must also be smaller than cell number.
  # 邻居数量不能超过实际cell数量。

  K_NEIGHBORS_ACTUAL <- min(
    as.integer(K_NEIGHBORS),
    n_cells_available - 1L
  )


  UMAP_N_NEIGHBORS_ACTUAL <- min(
    as.integer(UMAP_N_NEIGHBORS),
    n_cells_available - 1L
  )


  K_NEIGHBORS_ACTUAL <- max(
    2L,
    K_NEIGHBORS_ACTUAL
  )


  UMAP_N_NEIGHBORS_ACTUAL <- max(
    2L,
    UMAP_N_NEIGHBORS_ACTUAL
  )


  cat(
    "\n============================================================\n",
    "M10 PCA / CLUSTERING / UMAP\n",
    "============================================================\n",
    "Analysis assay: ", ANALYSIS_ASSAY, "\n",
    "Cells: ", n_cells_available, "\n",
    "HVGs: ", n_hvg_available, "\n",
    "PCs computed: ", N_PCS_ACTUAL, "\n",
    "PCs used: ", N_PCS_USE_ACTUAL, "\n",
    "KNN k: ", K_NEIGHBORS_ACTUAL, "\n",
    "UMAP neighbors: ", UMAP_N_NEIGHBORS_ACTUAL, "\n",
    "============================================================\n",
    sep = ""
  )


  # ==========================================================
  # 4. PCA
  # 4. 主成分分析
  # ==========================================================

  set.seed(
    RANDOM_SEED
  )


  seurat_clustered <- RunPCA(
    seurat_norm,
    assay = ANALYSIS_ASSAY,
    features = hvg_use,
    npcs = N_PCS_ACTUAL,
    seed.use = RANDOM_SEED,
    verbose = TRUE
  )


  # ==========================================================
  # 5. ELBOW PLOT
  # 5. PCA elbow图
  # ==========================================================

  p_elbow <- ElbowPlot(
    seurat_clustered,
    ndims = N_PCS_ACTUAL
  )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "PCA_elbow.png"
    ),
    p_elbow,
    width = 8,
    height = 6,
    dpi = 300
  )


  # ==========================================================
  # 6. NEIGHBOR GRAPH
  # 6. 构建邻居图
  # ==========================================================

  seurat_clustered <- FindNeighbors(
    seurat_clustered,
    reduction = "pca",
    dims = dims_use,
    k.param = K_NEIGHBORS_ACTUAL,

    # Explicit graph names make downstream clustering independent
    # of the assay name (RNA/SCT/etc.).
    #
    # 显式指定graph名称，使后续聚类不依赖RNA、SCT等assay名称。
    graph.name = c(
      CLUSTER_NN_GRAPH,
      CLUSTER_SNN_GRAPH
    ),

    verbose = TRUE
  )


  # ==========================================================
  # 7. GRAPH CLUSTERING AND RESOLUTION SENSITIVITY
  # 7. 图聚类与resolution敏感性检查
  # ==========================================================
  #
  # Multiple resolutions are evaluated only to inspect whether
  # cluster structure is stable across reasonable parameter choices.
  #
  # 多个resolution仅用于检查cluster划分对参数是否稳定。
  #
  # IMPORTANT:
  # clustree is a diagnostic tool. It does NOT automatically select
  # the biologically correct clustering resolution.
  #
  # 重要：
  # clustree只是诊断工具，不会自动替我们决定最佳resolution。
  # 最终正式分析仍明确使用CLUSTER_RESOLUTION。
  # ==========================================================

  resolution_grid <- sort(
    unique(
      c(
        CLUSTER_RESOLUTION_GRID,
        CLUSTER_RESOLUTION
      )
    )
  )


  # Run clustering on the explicitly named SNN graph.
  # 在前面明确命名的SNN graph上进行聚类。

  seurat_clustered <- FindClusters(
    seurat_clustered,
    graph.name = CLUSTER_SNN_GRAPH,
    resolution = resolution_grid,
    random.seed = RANDOM_SEED,
    verbose = TRUE
  )


  # Seurat stores clustering results using:
  #
  #   <graph.name>_res.<resolution>
  #
  # 因此根据graph名称自动生成resolution metadata前缀，
  # 而不是硬编码"RNA_snn_res."。

  resolution_prefix <- paste0(
    CLUSTER_SNN_GRAPH,
    "_res."
  )


  # ==========================================================
  # 7A. CLUSTERING-RESOLUTION DIAGNOSTIC
  # 7A. 聚类resolution稳定性诊断
  # ==========================================================
  #
  # clustree is used only as a diagnostic visualization.
  # It does NOT determine the final clustering resolution.
  #
  # clustree仅用于可视化不同resolution之间cluster的变化关系，
  # 不参与最终resolution的选择。
  #
  # IMPORTANT:
  # Some combinations of clustree, ggraph and ggplot2 may be
  # incompatible and can produce errors such as:
  #
  #   Unknown guide: edge_colourbar
  #
  # Such an error is a plotting/package-compatibility issue,
  # not a clustering failure. Therefore, clustree plotting is
  # wrapped in tryCatch() so that failure of this optional
  # diagnostic does not terminate the analysis pipeline.
  #
  # 注意：
  # 某些clustree / ggraph / ggplot2版本组合可能存在兼容性问题，
  # 例如出现：
  #
  #   Unknown guide: edge_colourbar
  #
  # 该错误仅代表辅助诊断图绘制失败，并不代表聚类失败。
  # 因此这里使用tryCatch()保护主分析流程。
  # ==========================================================

  clustree_success <- FALSE

  tryCatch(
    {

      p_clustree <- clustree::clustree(
        seurat_clustered[[]],
        prefix = resolution_prefix
      )

      p_clustree <- make_plot_background_transparent(
        p_clustree
      )

      clustree_file <- file.path(
        MDIR,
        "figures",
        "clustering_resolution_clustree.png"
      )


      # Explicitly build the ggraph/ggplot object and draw it through the ragg
      # device. This prevents a valid but completely blank PNG from being
      # produced by delayed grid rendering in some interactive RStudio sessions.
      #
      # 显式构建ggraph/ggplot图形对象，并通过ragg设备写入图片。这样可以避免
      # 部分RStudio交互式会话因grid延迟渲染而生成“文件存在但内容全白”的PNG。

      if (!requireNamespace("ragg", quietly = TRUE)) {
        stop(
          paste0(
            "Package 'ragg' is required to save the Clustree figure. ",
            "保存Clustree图片需要安装ragg包。"
          ),
          call. = FALSE
        )
      }


      ragg::agg_png(
        filename = clustree_file,
        width = 12,
        height = 8,
        units = "in",
        res = 300,
        background = "transparent"
      )


      tryCatch(
        {

          clustree_grob <- ggplot2::ggplotGrob(
            p_clustree
          )

          grid::grid.newpage()

          grid::grid.draw(
            clustree_grob
          )

        },
        finally = {

          if (grDevices::dev.cur() > 1L) {
            grDevices::dev.off()
          }

        }
      )

      clustree_success <- TRUE

      message(
        "Clustree resolution diagnostic saved successfully."
      )

    },
    error = function(e) {

      warning(
        paste0(
          "Clustree resolution diagnostic could not be generated. ",
          "This does NOT affect clustering results or downstream analysis.\n",
          "Reason: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )

    }
  )


  # ==========================================================
  # RESET FINAL CLUSTER IDENTITY
  # 设置正式cluster identity
  # ==========================================================
  #
  # Sensitivity analysis created several clustering columns.
  # Only CLUSTER_RESOLUTION is used as the official clustering
  # for downstream annotation and analysis.
  #
  # 上面会产生多个resolution对应的cluster结果。
  # 下游M11-M15只使用CLUSTER_RESOLUTION指定的正式结果。
  # ==========================================================

  final_resolution_col <- paste0(
    resolution_prefix,
    CLUSTER_RESOLUTION
  )


  if (
    !final_resolution_col %in%
    colnames(seurat_clustered[[]])
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Cannot find clustering column for final resolution: ",
        final_resolution_col,
        "\nAvailable clustering columns: ",
        paste(
          grep(
            "_res\\.",
            colnames(seurat_clustered[[]]),
            value = TRUE
          ),
          collapse = ", "
        )
      )
    )
  }


  Idents(
    seurat_clustered
  ) <- seurat_clustered[[]][[
    final_resolution_col
  ]]


  # ==========================================================
  # 8. UMAP
  # 8. UMAP降维
  # ==========================================================

  seurat_clustered <- RunUMAP(
    seurat_clustered,
    reduction = "pca",
    dims = dims_use,
    n.neighbors = UMAP_N_NEIGHBORS_ACTUAL,
    min.dist = UMAP_MIN_DIST,
    seed.use = RANDOM_SEED,
    verbose = TRUE
  )


  # Store cluster identity explicitly in metadata and preserve a natural
  # numerical order for all downstream tables, legends and grouped plots.
  #
  # If cluster IDs are kept as character values, R uses lexicographic order:
  #   0, 1, 10, 11, ..., 2, 20, ...
  # This is visually confusing but does not represent a change in clustering.
  # Numeric cluster IDs are therefore converted to an ordered factor whose
  # levels are 0, 1, 2, ..., while non-numeric IDs retain alphabetical order.
  #
  # 将cluster ID明确写入metadata，并为所有下游表格、图例和分组图片保留自然数值
  # 顺序。若cluster以character保存，R会按字典顺序排列为：
  #   0, 1, 10, 11, ..., 2, 20, ...
  # 这种顺序只会造成显示混乱，并不表示聚类结果发生变化。因此，纯数字cluster会被
  # 转换为levels按0、1、2……排列的factor；非数字ID则保留字母顺序。

  cluster_values <- as.character(
    Idents(
      seurat_clustered
    )
  )


  if (all(grepl("^[0-9]+$", cluster_values))) {

    cluster_levels <- as.character(
      sort(
        unique(
          as.integer(cluster_values)
        )
      )
    )

  } else {

    cluster_levels <- sort(
      unique(cluster_values)
    )
  }


  seurat_clustered$cluster <- factor(
    cluster_values,
    levels = cluster_levels
  )


  Idents(
    seurat_clustered
  ) <- seurat_clustered$cluster


  # ==========================================================
  # 9. CHECK REQUIRED METADATA FOR VISUALIZATION
  # 9. 检查绘图metadata
  # ==========================================================

  required_m10_meta <- c(
    SAMPLE_ID_COL,
    BIOLOGICAL_REPLICATE_COL,
    CONDITION_COL,
    "cluster"
  )


  missing_m10_meta <- setdiff(
    required_m10_meta,
    colnames(
      seurat_clustered[[]]
    )
  )


  if (length(missing_m10_meta) > 0L) {

    safe_stop_module(
      MODULE,
      paste0(
        "Missing metadata in M10 / M10缺少metadata: ",
        paste(
          missing_m10_meta,
          collapse = ", "
        )
      )
    )
  }


  # ==========================================================
  # 10. UMAP DIAGNOSTIC PLOTS
  # 10. UMAP诊断图
  # ==========================================================

  p_cluster <- DimPlot(
    seurat_clustered,
    reduction = "umap",
    group.by = "cluster",
    label = TRUE,
    repel = TRUE
  ) +
    ggtitle(
      "UMAP by cluster"
    )


  p_sample <- DimPlot(
    seurat_clustered,
    reduction = "umap",
    group.by = SAMPLE_ID_COL
  ) +
    ggtitle(
      "UMAP by sample"
    )


  p_condition <- DimPlot(
    seurat_clustered,
    reduction = "umap",
    group.by = CONDITION_COL
  ) +
    ggtitle(
      "UMAP by condition"
    )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "UMAP_overview.png"
    ),
    p_cluster +
      p_sample +
      p_condition,
    width = 20,
    height = 7,
    dpi = 300
  )

  # ==========================================================
  # 10.1 UMAP SPLIT BY SAMPLE
  # 10.1 按sample拆分UMAP
  # ==========================================================
  #
  # This plot shows the same clustering structure separately
  # for each sample.
  #
  # 该图把相同的UMAP结构按sample分别显示，
  # 用于检查某些cluster是否只由单一样本贡献，
  # 以及是否存在明显的sample-specific structure。
  #
  # This is a diagnostic visualization only.
  # It is NOT evidence of condition-level significance.
  #
  # 该图仅用于诊断，不能作为condition差异的统计检验。
  # ==========================================================

  p_sample_split <- DimPlot(
    seurat_clustered,
    reduction = "umap",
    group.by = "cluster",
    split.by = SAMPLE_ID_COL,
    label = FALSE
  ) +
    ggtitle(
      "UMAP split by sample"
    )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "UMAP_split_by_sample.png"
    ),
    p_sample_split,
    width = 4 * length(
      unique(
        seurat_clustered[[]][[SAMPLE_ID_COL]]
      )
    ),
    height = 6,
    dpi = 300,
    limitsize = FALSE
  )


  # ==========================================================
  # 10.2 CELL-CYCLE UMAP
  # 10.2 细胞周期UMAP
  # ==========================================================
  #
  # Cell-cycle phase is visualized on the existing UMAP
  # to assess whether major clusters are strongly associated
  # with proliferative state.
  #
  # 在现有UMAP上显示细胞周期状态，
  # 用于判断主要cluster是否明显受到细胞周期驱动。
  #
  # Cell-cycle scores are diagnostic covariates here.
  # They are not automatically regressed out.
  #
  # 此处细胞周期仅作为诊断变量，
  # 不默认从表达矩阵中回归去除。
  # ==========================================================

  if (
    "Phase" %in%
    colnames(
      seurat_clustered[[]]
    )
  ) {

    p_cell_cycle <- DimPlot(
      seurat_clustered,
      reduction = "umap",
      group.by = "Phase"
    ) +
      ggtitle(
        "UMAP by cell-cycle phase"
      )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "UMAP_cell_cycle.png"
      ),
      p_cell_cycle,
      width = 8,
      height = 6,
      dpi = 300
    )
  }

  # ==========================================================
  # 11. CLUSTER × SAMPLE SUMMARY
  # 11. cluster × sample汇总
  # ==========================================================

  cluster_sample <- seurat_clustered[[]] |>
    tibble::rownames_to_column(
      "cell"
    ) |>
    dplyr::count(
      cluster,
      .data[[SAMPLE_ID_COL]],
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]],
      name = "n_cells"
    )


  data.table::fwrite(
    cluster_sample,
    file.path(
      MDIR,
      "tables",
      "cluster_by_sample.csv"
    )
  )


  # ==========================================================
  # 12. SAVE
  # 12. 保存
  # ==========================================================

  save_checkpoint(
    seurat_clustered,
    "M10_clustered"
  )


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M10B. OPTIONAL MULTI-SAMPLE INTEGRATION / 可选多样本整合
# ============================================================
#
# 【Module scope / 模块范围】
# This optional sensitivity branch constructs an integrated low-dimensional
# representation when technical sample effects obscure shared cell identities.
# The primary non-integrated object is preserved as the reference analysis.
# 当技术性样本差异掩盖共同细胞身份时，本可选敏感性分支建立整合后的低维表示；
# 主分析的未整合对象仍被保留作为参考。
#
# 【Quality gate / 质量门槛】
# Integration should improve alignment of biologically comparable populations
# without erasing condition-specific states or merging distinct lineages. Compare
# markers, sample mixing and cluster structure before and after integration.
# 整合应改善生物学可比群体之间的对齐，同时不能抹去condition特异状态或合并不同
# 谱系。应比较整合前后的marker、样本混合程度和cluster结构。
#
# 【Inference restriction / 推断限制】
# Integrated values are not raw counts and are excluded from M13. All formal DE
# models continue to use unintegrated RNA counts aggregated by replicate.
# 整合值不是原始计数，不进入M13；所有正式DE模型仍使用按重复聚合的未整合RNA计数。
#
# 【默认为什么关闭 / Why disabled by default】
#
# sample差异既可能来自技术batch，也可能来自真实疾病效应。
#
# Sample differences may represent either technical batch
# effects or genuine condition-associated biology.
#
# 因此不能仅仅为了让UMAP“混得更漂亮”而进行integration。
#
# Integration should NOT be performed merely to make samples
# visually mix better.
#
#
# 【原则 / Principle】
#
# integrated representation：
#
#   可用于cluster / visualization / annotation
#
# raw RNA counts：
#
#   用于M13 pseudobulk differential expression
#
# 正式condition-level DE绝不能使用integrated expression。
#
#
# 【什么时候打开 / When to enable】
#
# 只有M10显示明显技术batch effect，
# 并且有充分实验设计依据时才建议开启。
#
#
# 【Output / 输出】
#
# integrated_UMAP.png
# M10B_integrated_visualization_only.rds
#
# ============================================================


if (RUN_M10B_INTEGRATION) {

  MODULE <- "M10B_integration"
  MDIR <- module_dir(MODULE)


  seurat_integrated <- load_checkpoint(
    "M10_clustered"
  )


  if (
    INTEGRATION_METHOD != "CCA"
  ) {

    safe_stop_module(
      MODULE,
      "Current template only implements CCA integration."
    )
  }


  DefaultAssay(
    seurat_integrated
  ) <- "RNA"

  if (!inherits(seurat_integrated[["RNA"]], "Assay5")) {
    safe_stop_module(
      MODULE,
      paste0(
        "M10B requires a Seurat v5 Assay5 RNA assay. ",
        "M10B需要Seurat v5 Assay5格式的RNA assay。"
      )
    )
  }


  # ==========================================================
  # 1. PREPARE SEURAT v5 SAMPLE-SPECIFIC LAYERS
  # ==========================================================
  #
  # IntegrateLayers() expects sample-specific layers.
  #
  # M09 previously joined RNA layers for Bioconductor
  # compatibility. Therefore M10B must split the RNA assay
  # again by sample before integration.
  #
  # M09为了OSCA分析已经把RNA layers合并。
  # 这里进行Seurat v5 integration时，需要重新按sample split。
  # ==========================================================

  if (
    !SAMPLE_ID_COL %in%
    colnames(
      seurat_integrated[[]]
    )
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Missing sample column for integration: ",
        SAMPLE_ID_COL
      )
    )
  }


  if (
    inherits(
      seurat_integrated[["RNA"]],
      "Assay5"
    )
  ) {

    current_layers <- Layers(
      seurat_integrated[["RNA"]]
    )


    # Only split if RNA is currently joined.
    # 只有当前没有sample-specific counts layers时才split。

    sample_specific_counts <- grep(
      "^counts\\.",
      current_layers,
      value = TRUE
    )


    if (length(sample_specific_counts) <= 1L) {

      sample_factor <- factor(
        as.character(
          seurat_integrated[[]][[SAMPLE_ID_COL]]
        )
      )


      seurat_integrated[["RNA"]] <- split(
        seurat_integrated[["RNA"]],
        f = sample_factor
      )
    }
  }


  # ==========================================================
  # 2. RE-NORMALIZE SAMPLE-SPECIFIC LAYERS FOR INTEGRATION
  # ==========================================================
  #
  # This normalization is ONLY for the optional Seurat
  # integration representation.
  #
  # 这里重新NormalizeData并不是替代M09的OSCA normalization。
  #
  # 它只服务于M10B的CCA integration。
  # ==========================================================

  seurat_integrated <- NormalizeData(
    seurat_integrated,
    verbose = FALSE
  )


  seurat_integrated <- FindVariableFeatures(
    seurat_integrated,
    selection.method = "vst",
    nfeatures = min(
      as.integer(N_HVG),
      nrow(seurat_integrated)
    ),
    verbose = FALSE
  )


  seurat_integrated <- ScaleData(
    seurat_integrated,
    features = VariableFeatures(
      seurat_integrated
    ),
    verbose = FALSE
  )


  n_pcs_int <- min(
    as.integer(N_PCS_COMPUTE),
    length(
      VariableFeatures(
        seurat_integrated
      )
    ) - 1L,
    ncol(seurat_integrated) - 1L
  )


  if (n_pcs_int < 2L) {

    safe_stop_module(
      MODULE,
      "Too few cells/HVGs for integration PCA."
    )
  }


  seurat_integrated <- RunPCA(
    seurat_integrated,
    features = VariableFeatures(
      seurat_integrated
    ),
    npcs = n_pcs_int,
    seed.use = RANDOM_SEED,
    verbose = FALSE
  )


  # ==========================================================
  # 3. CCA INTEGRATION
  # ==========================================================

  seurat_integrated <- IntegrateLayers(
    object = seurat_integrated,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    verbose = TRUE
  )


  n_pcs_int_use <- min(
    as.integer(N_PCS_USE),
    n_pcs_int
  )


  dims_int <- seq_len(
    n_pcs_int_use
  )


  # ==========================================================
  # 4. GRAPH / CLUSTER / UMAP ON INTEGRATED SPACE
  # ==========================================================

  seurat_integrated <- FindNeighbors(
    seurat_integrated,
    reduction = "integrated.cca",
    dims = dims_int
  )


  seurat_integrated <- FindClusters(
    seurat_integrated,
    resolution = CLUSTER_RESOLUTION,
    random.seed = RANDOM_SEED
  )


  seurat_integrated <- RunUMAP(
    seurat_integrated,
    reduction = "integrated.cca",
    dims = dims_int,
    reduction.name = "umap.integrated",
    n.neighbors = min(
      as.integer(UMAP_N_NEIGHBORS),
      ncol(seurat_integrated) - 1L
    ),
    min.dist = UMAP_MIN_DIST,
    seed.use = RANDOM_SEED
  )


  seurat_integrated$integrated_cluster <-
    as.character(
      Idents(
        seurat_integrated
      )
    )


  # ==========================================================
  # 5. INTEGRATION DIAGNOSTIC
  # ==========================================================

  p_int_sample <- DimPlot(
    seurat_integrated,
    reduction = "umap.integrated",
    group.by = SAMPLE_ID_COL
  ) +
    ggtitle(
      "Integrated UMAP by sample"
    )


  p_int_condition <- DimPlot(
    seurat_integrated,
    reduction = "umap.integrated",
    group.by = CONDITION_COL
  ) +
    ggtitle(
      "Integrated UMAP by condition"
    )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "integrated_UMAP.png"
    ),
    p_int_sample +
      p_int_condition,
    width = 14,
    height = 6,
    dpi = 300
  )


  # ==========================================================
  # INTEGRATED UMAP SPLIT BY SAMPLE
  # 按sample拆分integrated UMAP
  # ==========================================================
  #
  # This plot provides a direct sample-by-sample view of the
  # integrated representation.
  #
  # 该图逐个sample显示integration后的UMAP，
  # 用于辅助判断integration后各sample是否仍保留
  # comparable cell populations and clustering structure。
  #
  # Integration is used for visualization/clustering only.
  # Formal condition-level DE still uses raw RNA counts.
  #
  # Integration仅用于visualization/clustering。
  # 正式condition-level DE仍使用raw RNA counts。
  # ==========================================================

  p_int_split <- DimPlot(
    seurat_integrated,
    reduction = "umap.integrated",
    group.by = "seurat_clusters",
    split.by = SAMPLE_ID_COL,
    label = FALSE
  ) +
    ggtitle(
      "Integrated UMAP split by sample"
    )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "integrated_UMAP_split_by_sample.png"
    ),
    p_int_split,
    width = 4 * length(
      unique(
        seurat_integrated[[]][[SAMPLE_ID_COL]]
      )
    ),
    height = 6,
    dpi = 300,
    limitsize = FALSE
  )

  save_checkpoint(
    seurat_integrated,
    "M10B_integrated_visualization_only"
  )


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M11. CLUSTER MARKERS AND CELL-TYPE ANNOTATION / Cluster marker与细胞类型注释
# ============================================================
#
# 【Module scope / 模块范围】
# Cluster-enriched genes are identified to characterize transcriptional identity,
# canonical liver marker panels are visualized, and reviewed cell-type labels are
# written back to the object through an explicit manual annotation table.
# 本模块识别cluster富集基因以描述转录身份，可视化经典肝脏marker组合，并通过明确的
# 人工注释表将审阅后的细胞类型标签写回对象。
#
# 【Quality gate / 质量门槛】
# Each label should be supported by multiple positive markers, absence of conflicting
# lineage markers, cluster-level expression patterns and liver context. Ambiguous,
# mixed or low-quality clusters should retain a cautious label rather than receive
# an over-specific assignment.
# 每个标签应由多个阳性marker、冲突谱系marker的缺失、cluster层面表达模式和肝脏背景
# 共同支持。对模糊、混合或低质量cluster，应保留谨慎标签，而不应过度细分。
#
# 【Inference restriction / 推断限制】
# FindAllMarkers compares cells between clusters for annotation and is exploratory.
# It does not test disease effects and must not be reported as a replicate-level
# case-control result.
# FindAllMarkers用于cluster间细胞比较和注释，属于探索性分析；它不检验疾病效应，
# 不能作为生物学重复层面的病例—对照结果报告。
#
# 【Purpose / 目的】
#
# 找出每个cluster的marker genes，
# 并结合经典肝脏marker进行人工cell-type annotation。
#
#
# 【非常重要的统计区别】
#
# FindAllMarkers：
#
#   cell-level cluster A vs other clusters
#   用于解释cluster identity和annotation
#
# M13 edgeR：
#
#   Case vs Control
#   biological-replicate-level pseudobulk
#
# 两者完全不是同一个问题。
#
#
# 【注释原则 / Annotation principles】
#
# 1. 不依赖单个marker。
# 2. 使用多个一致marker。
# 3. 同时检查冲突marker。
# 4. 结合top markers + DotPlot + UMAP。
# 5. 疾病可能改变经典marker表达。
# 6. cluster编号本身没有生物学意义。
#
#
# 【Output / 输出】
#
# all_cluster_markers.csv
# top10_markers_per_cluster.csv
# top_marker_heatmap.png
# top_marker_heatmap_grouped_by_annotation.png
# cell_type_marker_ComplexHeatmap.png
# ComplexHeatmap_marker_panel.csv
# liver_marker_dotplot.png
# canonical_marker_featureplot.png
# manual_cluster_annotation.csv
# UMAP_cell_types.png
# UMAP_cell_types_split_by_sample.png
# M11_annotated.rds
# ============================================================


if (RUN_M11_MARKERS_ANNOTATION) {

  MODULE <- "M11_markers_annotation"
  MDIR <- module_dir(MODULE)


  seurat_clustered <- load_checkpoint(
    "M10_clustered"
  )


  ANALYSIS_ASSAY <- attr(
    seurat_clustered,
    "ANALYSIS_ASSAY"
  )


  if (
    is.null(ANALYSIS_ASSAY) ||
    !ANALYSIS_ASSAY %in% SeuratObject::Assays(seurat_clustered)
  ) {

    ANALYSIS_ASSAY <- DefaultAssay(
      seurat_clustered
    )
  }


  DefaultAssay(
    seurat_clustered
  ) <- ANALYSIS_ASSAY


  # ==========================================================
  # 1. PREPARE ASSAY FOR MARKER TESTING
  # ==========================================================
  #
  # FindAllMarkers needs a unified expression layer.
  #
  # 对RNA Assay5，如果仍存在多个data layers，
  # 先join后再进行marker analysis。
  # ==========================================================

  if (
    ANALYSIS_ASSAY == "RNA" &&
    inherits(
      seurat_clustered[["RNA"]],
      "Assay5"
    )
  ) {

    seurat_clustered <- join_assay_layers_if_needed(
      seurat_clustered,
      assay = "RNA",
      verbose = TRUE
    )
  }


  if (
    !"cluster" %in%
    colnames(
      seurat_clustered[[]]
    )
  ) {

    safe_stop_module(
      MODULE,
      "cluster metadata is missing."
    )
  }


  # Re-establish natural cluster order inside M11 as well as M10. This makes
  # M11 safe to rerun from an older M10 checkpoint in which cluster metadata
  # may still be stored as character values.
  #
  # M11会再次建立自然cluster顺序，因此即使从cluster仍以character保存的旧M10
  # checkpoint重新运行，也能保证热图、DotPlot和图例按照数值顺序显示。

  cluster_values <- as.character(
    seurat_clustered$cluster
  )


  if (all(grepl("^[0-9]+$", cluster_values))) {

    cluster_levels <- as.character(
      sort(
        unique(
          as.integer(cluster_values)
        )
      )
    )

  } else {

    cluster_levels <- sort(
      unique(cluster_values)
    )
  }


  seurat_clustered$cluster <- factor(
    cluster_values,
    levels = cluster_levels
  )


  Idents(
    seurat_clustered
  ) <- seurat_clustered$cluster


  # ==========================================================
  # 2. FIND CLUSTER MARKERS
  # ==========================================================

  all_markers <- FindAllMarkers(
    seurat_clustered,
    assay = ANALYSIS_ASSAY,
    only.pos = TRUE,
    test.use = "wilcox",
    min.pct = 0.20,
    logfc.threshold = 0.25,
    return.thresh = 0.05,
    verbose = TRUE
  )


  if (
    is.null(all_markers) ||
    nrow(all_markers) == 0L
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "FindAllMarkers returned no markers. ",
        "FindAllMarkers没有得到marker。"
      )
    )
  }


  data.table::fwrite(
    all_markers,
    file.path(
      MDIR,
      "tables",
      "all_cluster_markers.csv"
    )
  )


  # ==========================================================
  # 3. TOP 10 MARKERS PER CLUSTER
  # ==========================================================

  logfc_candidates <- c(
    "avg_log2FC",
    "avg_logFC"
  )


  logfc_col <- logfc_candidates[
    logfc_candidates %in%
      colnames(all_markers)
  ][1]


  if (
    is.na(logfc_col)
  ) {

    safe_stop_module(
      MODULE,
      "Cannot identify log-fold-change column in FindAllMarkers output."
    )
  }


  top10 <- all_markers |>
    dplyr::group_by(
      cluster
    ) |>
    dplyr::slice_max(
      order_by = .data[[logfc_col]],
      n = 10,
      with_ties = FALSE
    ) |>
    dplyr::ungroup()


  data.table::fwrite(
    top10,
    file.path(
      MDIR,
      "tables",
      "top10_markers_per_cluster.csv"
    )
  )

  # ==========================================================
  # TOP-MARKER HEATMAP
  # Top marker热图
  # ==========================================================
  #
  # The heatmap visualizes representative cluster-enriched
  # genes identified by FindAllMarkers().
  #
  # 该热图显示FindAllMarkers()鉴定出的代表性cluster marker，
  # 用于检查不同cluster是否具有清晰且相互一致的表达特征。
  #
  # Marker heatmaps support annotation but do not replace
  # canonical-marker validation.
  #
  # Marker热图用于支持细胞注释，
  # 但不能替代经典marker的独立验证。
  # ==========================================================

  heatmap_genes <- unique(
    top10$gene
  )


  heatmap_genes <- intersect(
    heatmap_genes,
    rownames(
      seurat_clustered
    )
  )


  if (length(heatmap_genes) > 0L) {

    heatmap_object <- ScaleData(
      seurat_clustered,
      features = heatmap_genes,
      verbose = FALSE
    )


    p_marker_heatmap <- DoHeatmap(
      heatmap_object,
      features = heatmap_genes,
      group.by = "cluster"
    )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "top_marker_heatmap.png"
      ),
      p_marker_heatmap,
      width = 12,
      height = 10,
      dpi = 300
    )

    rm(heatmap_object)
  }


  # ==========================================================
  # 4. GENERIC LIVER MARKER PANEL
  # ==========================================================

  species_marker_key <- tolower(SPECIES)


  if (species_marker_key == "mouse") {

    marker_panel <- list(

      Hepatocyte = c(
        "Alb", "Ttr", "Apoa1", "Apoa2",
        "Asgr1", "Hnf4a"
      ),

      Cholangiocyte = c(
        "Krt19", "Krt8", "Krt18",
        "Epcam", "Sox9"
      ),

      Endothelial = c(
        "Pecam1", "Kdr", "Vwf", "Emcn"
      ),

      LSEC = c(
        "Stab2", "Clec4g", "Lyve1", "Kdr"
      ),

      Kupffer = c(
        "Clec4f", "Timd4", "Marco",
        "Vsig4", "C1qa", "C1qb"
      ),

      Monocyte_macrophage = c(
        "Lyz2", "Ccr2", "Lgals3",
        "Trem2", "Spp1"
      ),

      Stellate_fibroblast = c(
        "Col1a1", "Col1a2", "Col3a1",
        "Dcn", "Lum", "Lrat"
      ),

      Activated_fibroblast = c(
        "Col1a1", "Acta2", "Postn", "Timp1"
      ),

      T_cell = c(
        "Cd3d", "Cd3e", "Trbc1", "Trbc2"
      ),

      NK_cell = c(
        "Nkg7", "Klrd1", "Prf1", "Gzmb"
      ),

      B_cell = c(
        "Cd79a", "Ms4a1", "Cd74"
      ),

      Neutrophil = c(
        "S100a8", "S100a9", "Ly6g", "Csf3r"
      ),

      Dendritic = c(
        "Flt3", "Xcr1", "Clec10a", "H2-Ab1"
      )
    )


  } else if (species_marker_key == "human") {

    marker_panel <- list(

      Hepatocyte = c(
        "ALB", "TTR", "APOA1",
        "ASGR1", "HNF4A"
      ),

      Cholangiocyte = c(
        "KRT19", "KRT8", "KRT18",
        "EPCAM", "SOX9"
      ),

      Endothelial = c(
        "PECAM1", "KDR", "VWF", "EMCN"
      ),

      LSEC = c(
        "STAB2", "CLEC4G", "LYVE1"
      ),

      Kupffer_macrophage = c(
        "MARCO", "VSIG4", "C1QA", "C1QB"
      ),

      Monocyte_macrophage = c(
        "LYZ", "CCR2", "LGALS3",
        "TREM2", "SPP1"
      ),

      Stellate_fibroblast = c(
        "COL1A1", "COL1A2", "COL3A1",
        "DCN", "LUM", "LRAT"
      ),

      T_cell = c(
        "CD3D", "CD3E", "TRBC1"
      ),

      NK_cell = c(
        "NKG7", "KLRD1", "PRF1"
      ),

      B_cell = c(
        "CD79A", "MS4A1", "CD74"
      )
    )


  } else {

    safe_stop_module(
      MODULE,
      paste0(
        "A validated liver annotation marker panel is not yet ",
        "configured for SPECIES = '",
        SPECIES,
        "'.\n",
        "Please provide a species-specific marker panel before ",
        "manual cell-type annotation."
      )
    )
  }


  # ------------------------------------------------------------
  # Prepare liver marker panel for DotPlot
  # 为DotPlot准备肝脏marker panel
  # ------------------------------------------------------------
  #
  # Step 1:
  # Keep only genes that are actually present in the dataset.
  #
  # 第一步：
  # 只保留当前数据集中真实存在的marker genes。
  #
  # Some marker genes may legitimately occur in more than one
  # biological category. For example:
  #
  #   Kdr     -> Endothelial / LSEC
  #   Col1a1  -> Stellate_fibroblast / Activated_fibroblast
  #
  # 同一个marker出现在多个相关细胞类型中在生物学上是合理的，
  # 因此这里不修改原始marker_panel的生物学定义。
  # ------------------------------------------------------------

  marker_panel_present <- lapply(
    marker_panel,
    function(x) {

      x[
        x %in%
          rownames(
            seurat_clustered
          )
      ]

    }
  )


  # Remove empty marker categories.
  # 删除当前数据中一个marker都不存在的类别。

  marker_panel_present <-
    marker_panel_present[
      lengths(
        marker_panel_present
      ) > 0L
    ]


  # ------------------------------------------------------------
  # Remove duplicated genes specifically for DotPlot
  # 仅为DotPlot自动处理重复marker
  # ------------------------------------------------------------
  #
  # Seurat::DotPlot() internally converts feature names to a factor.
  # Factor levels must be unique.
  #
  # Seurat::DotPlot()内部会将feature名称转换成factor，
  # 而factor levels不能重复。
  #
  # Therefore, if the same marker occurs in multiple biological
  # categories, passing the original list directly to DotPlot()
  # may produce:
  #
  #   factor level [...] is duplicated
  #
  # 为保持母版的通用性：
  #
  #   1. 原始marker_panel保持不变；
  #   2. 只在用于DotPlot时生成一个去重后的版本；
  #   3. marker第一次出现时保留；
  #   4. 后续重复出现自动删除。
  #
  # This preserves biological information in marker_panel while
  # ensuring compatibility with Seurat DotPlot.
  # ------------------------------------------------------------

  marker_panel_dotplot <- list()

  genes_already_used <- character(0)


  for (
    marker_group in names(
      marker_panel_present
    )
  ) {

    current_genes <-
      marker_panel_present[[marker_group]]


    current_genes <- current_genes[
      !current_genes %in%
        genes_already_used
    ]


    if (
      length(
        current_genes
      ) > 0L
    ) {

      marker_panel_dotplot[[marker_group]] <- current_genes


      genes_already_used <- c(
        genes_already_used,
        current_genes
      )
    }
  }


  # Final safety check:
  # DotPlot features must be globally unique.
  #
  # 最终安全检查：
  # 确保用于DotPlot的所有marker在全局范围内唯一。

  dotplot_features <-
    unlist(
      marker_panel_dotplot,
      use.names = FALSE
    )


  if (
    anyDuplicated(
      dotplot_features
    ) > 0L
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Duplicated marker genes remain after DotPlot preprocessing: ",
        paste(
          unique(
            dotplot_features[
              duplicated(
                dotplot_features
              )
            ]
          ),
          collapse = ", "
        )
      )
    )
  }


  # ------------------------------------------------------------
  # Draw marker DotPlot
  # 绘制marker DotPlot
  # ------------------------------------------------------------

  if (
    length(
      marker_panel_dotplot
    ) > 0L
  ) {

    p_dot <- DotPlot(
      seurat_clustered,
      features =
        marker_panel_dotplot,
      group.by = "cluster"
    ) +
      RotatedAxis() +
      # ggplot2 draws the first discrete y level at the bottom. Reversing the
      # display limits places cluster 0 at the top and the largest cluster ID
      # at the bottom, so the figure reads from 0 to 24 from top to bottom.
      # ggplot2默认把第一个离散y轴level画在底部。反转显示limits后，cluster 0位于
      # 顶部、最大cluster ID位于底部，因此图片从上到下按照0至24阅读。
      scale_y_discrete(
        limits = rev(cluster_levels),
        drop = FALSE
      )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "liver_marker_dotplot.png"
      ),
      p_dot,
      width = 20,
      height = 10,
      dpi = 300,
      limitsize = FALSE
    )
  }

  # ==========================================================
  # CANONICAL-MARKER FEATURE PLOTS
  # 经典marker FeaturePlot
  # ==========================================================
  #
  # Select one representative marker from each marker category.
  # Because marker_panel is species-specific, this automatically
  # produces species-compatible FeaturePlot genes.
  #
  # 每个细胞类型自动选择一个代表性marker。
  # 因为marker_panel本身已经按物种定义，
  # 因此这里不再硬编码mouse gene symbols。
  # ==========================================================

  featureplot_markers <- vapply(
    marker_panel_present,
    function(x) {

      if (length(x) == 0L) {
        return(NA_character_)
      }

      x[[1L]]
    },
    FUN.VALUE = character(1)
  )


  featureplot_markers <- unique(
    stats::na.omit(
      featureplot_markers
    )
  )


  featureplot_markers <- intersect(
    featureplot_markers,
    rownames(seurat_clustered)
  )


  if (length(featureplot_markers) > 0L) {

    p_feature <- FeaturePlot(
      seurat_clustered,
      features = featureplot_markers,
      reduction = "umap",
      ncol = min(
        4L,
        length(featureplot_markers)
      )
    )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "canonical_marker_featureplot.png"
      ),
      p_feature,
      width = 14,
      height = max(
        6,
        ceiling(
          length(featureplot_markers) / 4
        ) * 4
      ),
      dpi = 300,
      limitsize = FALSE
    )
  }



  # ==========================================================
  # 5. MANUAL ANNOTATION TEMPLATE
  # ==========================================================

  annotation_file <- file.path(
    MDIR,
    "tables",
    "manual_cluster_annotation.csv"
  )


  cluster_ids <- unique(
    as.character(
      seurat_clustered$cluster
    )
  )

  if (all(grepl("^[0-9]+$", cluster_ids))) {
    cluster_ids <- cluster_ids[
      order(
        as.integer(cluster_ids)
      )
    ]
  } else {
    cluster_ids <- sort(cluster_ids)
  }


  if (
    CREATE_MANUAL_ANNOTATION_TEMPLATE &&
    !file.exists(annotation_file)
  ) {

    annotation_template <- data.frame(
      cluster = cluster_ids,
      cell_type = "",
      confidence = "",
      evidence = "",
      comments = "",
      stringsAsFactors = FALSE
    )


    data.table::fwrite(
      annotation_template,
      annotation_file
    )


    pending_message <- paste0(
      "Manual annotation template generated / 已生成注释模板:\n",
      annotation_file,
      "\nPlease complete every cell_type and rerun M11. ",
      "请填写全部cell_type后重新运行M11。"
    )

    write_module_status(
      MODULE,
      "WAITING_FOR_ANNOTATION",
      pending_message
    )

    stop(
      pending_message,
      call. = FALSE
    )
  }


  # ==========================================================
  # 6. READ MANUAL ANNOTATION
  # ==========================================================

  if (
    file.exists(annotation_file)
  ) {

    anno <- data.table::fread(
      annotation_file,
      data.table = FALSE
    )


    if (
      !all(
        c(
          "cluster",
          "cell_type"
        ) %in%
        colnames(anno)
      )
    ) {

      safe_stop_module(
        MODULE,
        paste0(
          "manual_cluster_annotation.csv must contain ",
          "cluster and cell_type columns."
        )
      )
    }


    anno$cluster <- trimws(
      as.character(
        anno$cluster
      )
    )

    anno$cell_type <- trimws(
      as.character(
        anno$cell_type
      )
    )

    if (anyDuplicated(anno$cluster) > 0L) {
      safe_stop_module(
        MODULE,
        paste0(
          "manual_cluster_annotation.csv contains duplicated cluster IDs / ",
          "注释表包含重复cluster。"
        )
      )
    }

    missing_clusters <- setdiff(
      cluster_ids,
      anno$cluster
    )

    extra_clusters <- setdiff(
      anno$cluster,
      cluster_ids
    )

    if (
      length(missing_clusters) > 0L ||
      length(extra_clusters) > 0L
    ) {
      safe_stop_module(
        MODULE,
        paste0(
          "Annotation clusters do not match the current M10 object. ",
          "注释表与当前M10 cluster不一致。\n",
          "Missing / 缺少: ",
          paste(missing_clusters, collapse = ", "),
          "\nExtra / 多出: ",
          paste(extra_clusters, collapse = ", ")
        )
      )
    }

    if (any(is.na(anno$cell_type) | !nzchar(anno$cell_type))) {
      safe_stop_module(
        MODULE,
        paste0(
          "Every cluster must have a non-empty cell_type / ",
          "每个cluster都必须填写非空cell_type。"
        )
      )
    }


    map <- setNames(
      as.character(
        anno$cell_type
      ),
      as.character(
        anno$cluster
      )
    )


    assigned <- unname(
      map[
        as.character(
          seurat_clustered$cluster
        )
      ]
    )


    seurat_clustered$cell_type <-
      assigned
  }


  # Complete annotation is required before downstream modules.
  # 进入下游模块前必须完成全部人工注释。

  if (
    !"cell_type" %in%
    colnames(
      seurat_clustered[[]]
    )
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Complete manual annotation is required before M11 can finish / ",
        "完成人工注释后M11才能结束。"
      )
    )
  }


  # ==========================================================
  # 7. ANNOTATION-GROUPED TOP-20 MARKER HEATMAP
  # 7. 按细胞类型注释分组的Top 20 marker热图
  # ==========================================================
  #
  # The original top-marker heatmap remains grouped by numerical cluster
  # identity and retains the complete marker set. This additional summary
  # heatmap displays exactly 20 representative marker genes in 10 presentation
  # categories. The nine most abundant reviewed cell types are retained, while
  # all remaining cell types are combined into a tenth category named Others.
  # One leading marker is first retained for every displayed category; the
  # remaining slots are filled by the largest marker log-fold changes.
  #
  # 原始Top marker热图仍按照数值cluster identity展示，并保留完整marker集合。
  # 本新增汇总热图严格显示20个代表性marker genes和10个展示分类。保留cell数量
  # 最多的9种已审阅cell type，其余cell type统一合并为第10类Others。先为每个
  # 展示分类保留1个最高marker，再按照marker log-fold change由高到低补足剩余名额。
  #
  # The nine retained cell types are ordered by decreasing cell count, with
  # Others fixed as the tenth and final category. Cluster numbers are omitted.
  #
  # 保留的9种cell type按照cell数量从多到少排序，Others固定为第10个且位于最后。
  # 本图不再显示cluster编号，以提供更简洁的细胞类型层面展示。
  #
  # This grouping is used only inside the figure. It does not overwrite the
  # original cell_type metadata, merge clusters, change cell identities or alter
  # marker-test results.
  #
  # 该合并仅用于本图展示，不会覆盖原始cell_type metadata，不会合并cluster、改变
  # cell identity或修改marker检验结果。
  # ==========================================================

  if (length(heatmap_genes) > 0L) {


    cluster_annotation_order <- unique(
      data.frame(
        cluster = as.character(
          seurat_clustered$cluster
        ),
        cell_type = as.character(
          seurat_clustered$cell_type
        ),
        stringsAsFactors = FALSE
      )
    )


    cluster_annotation_order$cluster_factor <- factor(
      cluster_annotation_order$cluster,
      levels = cluster_levels
    )


    cluster_annotation_order <-
      cluster_annotation_order[
        order(
          cluster_annotation_order$cluster_factor
        ),
        ,
        drop = FALSE
      ]


    cell_type_counts <- sort(
      table(
        as.character(
          seurat_clustered$cell_type
        )
      ),
      decreasing = TRUE
    )


    major_cell_types <- names(
      utils::head(
        cell_type_counts,
        9L
      )
    )


    cluster_annotation_order$heatmap_cell_type <- ifelse(
      cluster_annotation_order$cell_type %in% major_cell_types,
      cluster_annotation_order$cell_type,
      "Others"
    )


    heatmap_cell_type_levels <- c(
      major_cell_types,
      "Others"
    )


    cluster_to_heatmap_cell_type <- setNames(
      cluster_annotation_order$heatmap_cell_type,
      cluster_annotation_order$cluster
    )


    annotation_marker_candidates <- top10 |>
      dplyr::mutate(
        cluster = as.character(cluster),
        heatmap_cell_type = unname(
          cluster_to_heatmap_cell_type[cluster]
        )
      ) |>
      dplyr::filter(
        !is.na(heatmap_cell_type),
        gene %in% rownames(seurat_clustered)
      )


    primary_annotation_markers <- annotation_marker_candidates |>
      dplyr::group_by(heatmap_cell_type) |>
      dplyr::slice_max(
        order_by = .data[[logfc_col]],
        n = 1L,
        with_ties = FALSE
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        heatmap_cell_type_order = factor(
          heatmap_cell_type,
          levels = heatmap_cell_type_levels
        )
      ) |>
      dplyr::arrange(heatmap_cell_type_order) |>
      dplyr::distinct(
        gene,
        .keep_all = TRUE
      ) |>
      dplyr::slice_head(
        n = 20L
      )


    n_additional_markers <- max(
      0L,
      20L - nrow(primary_annotation_markers)
    )


    additional_annotation_markers <- annotation_marker_candidates |>
      dplyr::filter(
        !gene %in% primary_annotation_markers$gene
      ) |>
      dplyr::arrange(
        dplyr::desc(
          .data[[logfc_col]]
        )
      ) |>
      dplyr::distinct(
        gene,
        .keep_all = TRUE
      ) |>
      dplyr::slice_head(
        n = n_additional_markers
      )


    annotation_heatmap_genes <- c(
      primary_annotation_markers$gene,
      additional_annotation_markers$gene
    )


    heatmap_cell_type_colours <- setNames(
      scales::hue_pal()(
        length(heatmap_cell_type_levels)
      ),
      heatmap_cell_type_levels
    )


    annotated_heatmap_object <- ScaleData(
      seurat_clustered,
      features = annotation_heatmap_genes,
      verbose = FALSE
    )


    annotated_heatmap_object$heatmap_cell_type <- factor(
      ifelse(
        as.character(
          annotated_heatmap_object$cell_type
        ) %in% major_cell_types,
        as.character(
          annotated_heatmap_object$cell_type
        ),
        "Others"
      ),
      levels = heatmap_cell_type_levels
    )


    p_marker_heatmap_grouped <- suppressMessages(
      DoHeatmap(
        annotated_heatmap_object,
        features = annotation_heatmap_genes,
        group.by = "heatmap_cell_type",
        group.colors = heatmap_cell_type_colours,
        label = TRUE,
        angle = 45,
        raster = TRUE
      ) +
        scale_fill_gradientn(
          colors = c(
            "white",
            "grey",
            "firebrick3"
          ),
          name = "Scaled\nexpression"
        ) +
        ggtitle(
          "Top 20 marker genes across 10 annotated categories"
        )
    )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "top_marker_heatmap_grouped_by_annotation.png"
      ),
      p_marker_heatmap_grouped,
      width = 14,
      height = 8,
      dpi = 300,
      limitsize = FALSE
    )


    # ========================================================
    # COMPLEXHEATMAP WITH CURATED MARKER BLOCKS
    # 使用人工策划marker分块的ComplexHeatmap
    # ========================================================
    #
    # This figure follows the matrix-first strategy used by pheatmap and
    # ComplexHeatmap rather than relying on the fixed layout of DoHeatmap.
    # Rows are split by the cell type that each marker is intended to support;
    # columns are split into the same nine major categories plus Others used in
    # the compact annotation heatmap above.
    #
    # 本图采用pheatmap/ComplexHeatmap的“先构建表达矩阵、再独立控制注释和分块”
    # 思路，不再受DoHeatmap固定布局限制。行按照marker所支持的细胞类型分块，列则
    # 沿用上方简化热图中的9个主要分类加Others，共10个展示分类。
    #
    # The normalized data layer generated by the selected M09 method is used.
    # Raw counts are deliberately excluded because direct log10(counts + 1)
    # would reintroduce library-size differences that OSCA/scran normalization
    # has already addressed. Expression is converted to a per-gene Z-score and
    # clipped symmetrically for visualization only.
    #
    # 本图使用M09所选标准化方法生成的normalized data layer。这里不使用raw counts，
    # 因为直接log10(counts + 1)会重新引入OSCA/scran已经处理的文库大小差异。表达值
    # 仅为绘图转换为逐基因Z-score，并进行对称截断。
    #
    # A fixed maximum number of cells is sampled from each displayed category.
    # This controls figure width and prevents abundant populations from visually
    # dominating the heatmap. The sampling is reproducible through RANDOM_SEED
    # and is used for visualization only.
    #
    # 每个展示分类最多抽取固定数量的cell，以控制图片宽度并避免高丰度群体在视觉上
    # 主导热图。抽样由RANDOM_SEED保证可重复，并且只影响绘图，不影响任何分析结果。
    # ========================================================

    if (
      !requireNamespace(
        "ComplexHeatmap",
        quietly = TRUE
      ) ||
      !requireNamespace(
        "circlize",
        quietly = TRUE
      )
    ) {
      safe_stop_module(
        MODULE,
        paste0(
          "ComplexHeatmap and circlize are required for the curated heatmap. ",
          "策划marker热图需要安装ComplexHeatmap和circlize。"
        )
      )
    }


    if (tolower(SPECIES) == "mouse") {

      complex_heatmap_marker_reference <- list(

        `B cell` = c(
          "Cd79a", "Cd79b", "Ms4a1", "Cd74", "Cd19"
        ),

        `Fibroblast / stellate cell` = c(
          "Col1a1", "Col1a2", "Col3a1", "Dcn", "Lum", "Lrat"
        ),

        `T cell` = c(
          "Cd3d", "Cd3e", "Cd3g", "Trac", "Trbc1"
        ),

        Monocyte = c(
          "Lyz2", "Ly6c2", "Ccr2", "Cd14", "Lgals3", "Chil3"
        ),

        pDC = c(
          "Siglech", "Bst2", "Ccr9", "Gzmb"
        ),

        `NK cell` = c(
          "Nkg7", "Klrd1", "Ncr1", "Prf1", "Xcl1"
        ),

        `LSEC / Endo L` = c(
          "Clec4g", "Stab2", "Lyve1", "Kdr", "Rspo3", "Oit3"
        ),

        Neutrophil = c(
          "S100a8", "S100a9", "Ly6g", "Csf3r", "Retnlg"
        ),

        cDC2 = c(
          "Clec10a", "Flt3", "Cd209a", "H2-Ab1", "Clec4b1"
        ),

        Cholangiocyte = c(
          "Krt19", "Krt8", "Krt18", "Epcam", "Sox9", "Mmp7"
        ),

        cDC1 = c(
          "Xcr1", "Clec9a", "Batf3", "Itgae", "Tlr3"
        ),

        Macrophage = c(
          "C1qa", "C1qb", "C1qc", "Adgre1", "Apoe", "Trem2"
        ),

        `Kupffer cell` = c(
          "Clec4f", "Timd4", "Marco", "Vsig4", "Cd5l", "Folr2"
        ),

        `Gamma-delta T cell` = c(
          "Trdc", "Tcrg-C1", "Rorc", "Cxcr6"
        ),

        Hepatocyte = c(
          "Alb", "Ttr", "Apoa1", "Apoa2", "Asgr1", "Hnf4a"
        )
      )

    } else if (tolower(SPECIES) == "human") {

      complex_heatmap_marker_reference <- list(

        `B cell` = c(
          "CD79A", "CD79B", "MS4A1", "CD74", "CD19"
        ),

        `Fibroblast / stellate cell` = c(
          "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "LRAT"
        ),

        `T cell` = c(
          "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1"
        ),

        Monocyte = c(
          "LYZ", "FCN1", "CCR2", "CD14", "LGALS3", "S100A8"
        ),

        pDC = c(
          "GZMB", "CLEC4C", "IL3RA", "TCF4"
        ),

        `NK cell` = c(
          "NKG7", "KLRD1", "NCR1", "PRF1", "XCL1"
        ),

        `LSEC / Endo L` = c(
          "CLEC4G", "STAB2", "LYVE1", "KDR", "RSPO3", "OIT3"
        ),

        Neutrophil = c(
          "S100A8", "S100A9", "CSF3R", "FCGR3B", "CXCR2"
        ),

        cDC2 = c(
          "CD1C", "CLEC10A", "FCER1A", "HLA-DRA", "FLT3"
        ),

        Cholangiocyte = c(
          "KRT19", "KRT8", "KRT18", "EPCAM", "SOX9", "MMP7"
        ),

        cDC1 = c(
          "XCR1", "CLEC9A", "BATF3", "ITGAE", "TLR3"
        ),

        Macrophage = c(
          "C1QA", "C1QB", "C1QC", "ADGRE1", "APOE", "TREM2"
        ),

        `Kupffer cell` = c(
          "TIMD4", "MARCO", "VSIG4", "CD5L", "FOLR2", "CD163"
        ),

        `Gamma-delta T cell` = c(
          "TRDC", "TRGC1", "RORC", "CXCR6"
        ),

        Hepatocyte = c(
          "ALB", "TTR", "APOA1", "APOA2", "ASGR1", "HNF4A"
        )
      )

    } else {

      complex_heatmap_marker_reference <- list()
    }


    all_reviewed_cell_types <- unique(
      cluster_annotation_order$cell_type
    )


    minor_cell_types <- setdiff(
      all_reviewed_cell_types,
      major_cell_types
    )


    complex_marker_rows <- list()
    complex_marker_index <- 1L


    for (
      display_category in heatmap_cell_type_levels
    ) {

      source_cell_types <- if (
        identical(
          display_category,
          "Others"
        )
      ) {
        minor_cell_types
      } else {
        display_category
      }


      for (source_cell_type in source_cell_types) {

        marker_values <-
          complex_heatmap_marker_reference[[
            source_cell_type
          ]]


        if (
          is.null(marker_values) ||
          length(marker_values) == 0L
        ) {

          source_clusters <-
            cluster_annotation_order$cluster[
              cluster_annotation_order$cell_type ==
                source_cell_type
            ]


          marker_values <- top10 |>
            dplyr::filter(
              as.character(cluster) %in% source_clusters
            ) |>
            dplyr::arrange(
              dplyr::desc(
                .data[[logfc_col]]
              )
            ) |>
            dplyr::distinct(
              gene,
              .keep_all = TRUE
            ) |>
            dplyr::slice_head(
              n = 5L
            ) |>
            dplyr::pull(gene)
        }


        complex_marker_rows[[complex_marker_index]] <-
          data.frame(
            marker_group = display_category,
            source_cell_type = source_cell_type,
            gene = marker_values,
            stringsAsFactors = FALSE
          )


        complex_marker_index <-
          complex_marker_index + 1L
      }
    }


    complex_marker_panel <- dplyr::bind_rows(
      complex_marker_rows
    ) |>
      dplyr::distinct(
        gene,
        .keep_all = TRUE
      ) |>
      dplyr::mutate(
        present = gene %in% rownames(
          seurat_clustered
        )
      )


    data.table::fwrite(
      complex_marker_panel,
      file.path(
        MDIR,
        "tables",
        "ComplexHeatmap_marker_panel.csv"
      )
    )


    complex_marker_panel_present <- complex_marker_panel |>
      dplyr::filter(present) |>
      dplyr::mutate(
        marker_group = factor(
          marker_group,
          levels = heatmap_cell_type_levels
        )
      ) |>
      dplyr::arrange(marker_group)


    if (nrow(complex_marker_panel_present) < 2L) {
      safe_stop_module(
        MODULE,
        paste0(
          "Too few curated markers were found for ComplexHeatmap. ",
          "ComplexHeatmap可用的策划marker过少。"
        )
      )
    }


    complex_heatmap_category <- factor(
      ifelse(
        as.character(
          seurat_clustered$cell_type
        ) %in% major_cell_types,
        as.character(
          seurat_clustered$cell_type
        ),
        "Others"
      ),
      levels = heatmap_cell_type_levels
    )


    names(complex_heatmap_category) <- colnames(
      seurat_clustered
    )


    set.seed(RANDOM_SEED)


    complex_cells_by_category <- lapply(
      heatmap_cell_type_levels,
      function(category_name) {

        category_cells <- names(
          complex_heatmap_category
        )[
          complex_heatmap_category == category_name
        ]


        if (
          length(category_cells) >
          COMPLEX_HEATMAP_MAX_CELLS_PER_CATEGORY
        ) {
          category_cells <- sample(
            category_cells,
            COMPLEX_HEATMAP_MAX_CELLS_PER_CATEGORY,
            replace = FALSE
          )
        }


        sort(category_cells)
      }
    )


    names(complex_cells_by_category) <-
      heatmap_cell_type_levels


    complex_heatmap_cells <- unlist(
      complex_cells_by_category,
      use.names = FALSE
    )


    if (length(complex_heatmap_cells) == 0L) {
      safe_stop_module(
        MODULE,
        paste0(
          "No cells were available for ComplexHeatmap. ",
          "ComplexHeatmap没有可用cell。"
        )
      )
    }


    complex_expression <- tryCatch(
      SeuratObject::LayerData(
        seurat_clustered,
        assay = ANALYSIS_ASSAY,
        layer = "data"
      ),
      error = function(e) {
        safe_stop_module(
          MODULE,
          paste0(
            "Unable to read normalized data for ComplexHeatmap / ",
            "无法读取ComplexHeatmap所需normalized data。\n",
            conditionMessage(e)
          )
        )
      }
    )


    complex_expression <- as.matrix(
      complex_expression[
        complex_marker_panel_present$gene,
        complex_heatmap_cells,
        drop = FALSE
      ]
    )


    complex_gene_means <- rowMeans(
      complex_expression
    )


    complex_gene_sds <- apply(
      complex_expression,
      1L,
      stats::sd
    )


    complex_variable_rows <-
      is.finite(complex_gene_sds) &
      complex_gene_sds > 0


    complex_expression <- complex_expression[
      complex_variable_rows,
      ,
      drop = FALSE
    ]


    complex_gene_means <- complex_gene_means[
      complex_variable_rows
    ]


    complex_gene_sds <- complex_gene_sds[
      complex_variable_rows
    ]


    complex_marker_panel_present <-
      complex_marker_panel_present[
        complex_variable_rows,
        ,
        drop = FALSE
      ]


    if (nrow(complex_expression) < 2L) {
      safe_stop_module(
        MODULE,
        paste0(
          "Too few variable markers remained for ComplexHeatmap. ",
          "ComplexHeatmap中保留的可变marker过少。"
        )
      )
    }


    complex_expression_z <- sweep(
      complex_expression,
      1L,
      complex_gene_means,
      FUN = "-"
    )


    complex_expression_z <- sweep(
      complex_expression_z,
      1L,
      complex_gene_sds,
      FUN = "/"
    )


    complex_expression_z[
      complex_expression_z >
        COMPLEX_HEATMAP_Z_LIMIT
    ] <- COMPLEX_HEATMAP_Z_LIMIT


    complex_expression_z[
      complex_expression_z <
        -COMPLEX_HEATMAP_Z_LIMIT
    ] <- -COMPLEX_HEATMAP_Z_LIMIT


    complex_column_split <- factor(
      as.character(
        complex_heatmap_category[
          complex_heatmap_cells
        ]
      ),
      levels = heatmap_cell_type_levels
    )


    complex_row_split <- factor(
      as.character(
        complex_marker_panel_present$marker_group
      ),
      levels = heatmap_cell_type_levels
    )


    complex_top_annotation <-
      ComplexHeatmap::HeatmapAnnotation(
        category = ComplexHeatmap::anno_block(
          gp = grid::gpar(
            fill = unname(
              heatmap_cell_type_colours[
                heatmap_cell_type_levels
              ]
            ),
            col = "white"
          ),
          labels = heatmap_cell_type_levels,
          labels_gp = grid::gpar(
            col = "black",
            fontsize = 9,
            fontface = "bold"
          )
        ),
        show_annotation_name = FALSE
      )


    complex_colour_function <- circlize::colorRamp2(
      c(
        -COMPLEX_HEATMAP_Z_LIMIT,
        0,
        COMPLEX_HEATMAP_Z_LIMIT
      ),
      c(
        "#2166AC",
        "white",
        "#B2182B"
      )
    )


    complex_heatmap_plot <- ComplexHeatmap::Heatmap(
      complex_expression_z,
      name = "Row Z-score",
      col = complex_colour_function,
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = FALSE,
      row_split = complex_row_split,
      column_split = complex_column_split,
      top_annotation = complex_top_annotation,
      column_title = paste0(
        "Curated marker expression across ",
        length(heatmap_cell_type_levels),
        " annotated categories"
      ),
      column_title_gp = grid::gpar(
        fontsize = 13,
        fontface = "bold"
      ),
      row_names_gp = grid::gpar(
        fontsize = 8
      ),
      row_title_gp = grid::gpar(
        fontsize = 9,
        fontface = "bold"
      ),
      row_gap = grid::unit(
        1.5,
        "mm"
      ),
      column_gap = grid::unit(
        1.5,
        "mm"
      ),
      border = TRUE,
      use_raster = TRUE,
      raster_quality = 2,
      heatmap_legend_param = list(
        title = "Row Z-score",
        at = c(
          -COMPLEX_HEATMAP_Z_LIMIT,
          0,
          COMPLEX_HEATMAP_Z_LIMIT
        ),
        legend_height = grid::unit(
          45,
          "mm"
        )
      )
    )


    complex_heatmap_file <- file.path(
      MDIR,
      "figures",
      "cell_type_marker_ComplexHeatmap.png"
    )


    ragg::agg_png(
      filename = complex_heatmap_file,
      width = 18,
      height = 14,
      units = "in",
      res = 300,
      background = "transparent"
    )


    tryCatch(
      {

        ComplexHeatmap::draw(
          complex_heatmap_plot,
          background = "transparent",
          heatmap_legend_side = "right",
          annotation_legend_side = "right"
        )

      },
      finally = {

        if (grDevices::dev.cur() > 1L) {
          grDevices::dev.off()
        }
      }
    )


    rm(annotated_heatmap_object)
  }


  # ==========================================================
  # 9. CELL-TYPE UMAP
  # ==========================================================

  p_celltype <- DimPlot(
    seurat_clustered,
    reduction = "umap",
    group.by = "cell_type",
    label = TRUE,
    repel = TRUE
  )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "UMAP_cell_types.png"
    ),
    p_celltype,
    width = 11,
    height = 8,
    dpi = 300
  )

  # ==========================================================
  # ANNOTATED UMAP SPLIT BY SAMPLE
  # 按sample拆分已注释UMAP
  # ==========================================================
  #
  # The final annotated cell types are displayed separately
  # for each sample.
  #
  # 将最终cell-type annotation按sample分别显示，
  # 用于检查主要细胞类型是否由多个样本共同支持，
  # 并识别潜在的sample-specific populations。
  #
  # This plot is particularly useful before composition and
  # pseudobulk analyses in M12-M13.
  #
  # 该图尤其适合在进入M12-M13的composition和
  # pseudobulk分析前进行最终检查。
  # ==========================================================

  p_cell_type_split <- DimPlot(
    seurat_clustered,
    reduction = "umap",
    group.by = "cell_type",
    split.by = SAMPLE_ID_COL,
    label = FALSE
  )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "UMAP_cell_types_split_by_sample.png"
    ),
    p_cell_type_split,
    width = 4 * length(
      unique(
        seurat_clustered[[]][[SAMPLE_ID_COL]]
      )
    ),
    height = 6,
    dpi = 300,
    limitsize = FALSE
  )

  validate_complete_annotation(
    seurat_clustered,
    MODULE
  )

  save_checkpoint(
    seurat_clustered,
    "M11_annotated"
  )


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M12. CELL-TYPE COMPOSITION / 细胞组成
# ============================================================
#
# 【Module scope / 模块范围】
# Cell counts and within-replicate proportions are summarized for every annotated
# cell type. The output shows representation patterns and supports data-quality
# review before any confirmatory abundance model is considered.
# 本模块汇总每个已注释细胞类型在各生物学重复中的细胞数和样本内比例，用于展示组成
# 模式，并在考虑验证性丰度模型前支持数据质量审查。
#
# 【Quality gate / 质量门槛】
# Confirm that all expected replicates are present and that apparent group shifts
# are not driven by one specimen, unequal cell recovery or upstream filtering.
# Zero counts are informative and should be distinguished from missing samples.
# 应确认所有预期重复均存在，并排除表面组间变化由单一样本、细胞捕获量不等或上游
# 过滤造成的可能。零计数具有信息，应与样本缺失明确区分。
#
# 【Inference restriction / 推断限制】
# These summaries are descriptive. Statistical claims about differential abundance
# require a dedicated replicate-aware compositional model and adequate replication.
# 本模块结果为描述性汇总。若要对差异丰度作统计结论，需要专门的重复感知组成模型
# 及足够的生物学重复。
#
# 【Purpose / 目的】
#
# 计算每个独立biological replicate中
# 不同cell type的cell数量和比例。
#
#
# 【统计单位 / Statistical unit】
#
# biological replicate，而不是单个cell。
#
# Biological replicate, NOT individual cells.
#
#
# 【为什么？ / Why?】
#
# 如果直接把一个condition中的所有cell合并，
# 捕获cell更多的动物会获得更大的统计权重，
# 从而产生pseudo-replication。
#
#
# 【注意 / Important】
#
# 当前M12只是descriptive composition。
#
# 它不是正式differential abundance test。
#
#
# 【Output / 输出】
#
# cell_type_composition_by_replicate.csv
# cell_type_composition_group_comparison_exploratory.csv
# cell_type_composition.png
# cell_type_composition_group_comparison_exploratory.png
#
# ============================================================


if (RUN_M12_COMPOSITION) {

  MODULE <- "M12_composition"
  MDIR <- module_dir(MODULE)


  obj <- load_checkpoint(
    "M11_annotated"
  )

  validate_complete_annotation(
    obj,
    MODULE
  )


  required_comp_meta <- c(
    BIOLOGICAL_REPLICATE_COL,
    CONDITION_COL,
    "cell_type"
  )


  missing_comp_meta <- setdiff(
    required_comp_meta,
    colnames(
      obj[[]]
    )
  )


  if (length(missing_comp_meta) > 0L) {

    safe_stop_module(
      MODULE,
      paste0(
        "Missing composition metadata: ",
        paste(
          missing_comp_meta,
          collapse = ", "
        )
      )
    )
  }

  replicate_condition_check <- obj[[]] |>
    dplyr::distinct(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]]
    ) |>
    dplyr::count(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      name = "n_conditions"
    )

  if (any(replicate_condition_check$n_conditions != 1L)) {
    safe_stop_module(
      MODULE,
      paste0(
        "A biological replicate maps to multiple conditions / ",
        "同一biological replicate对应多个condition。"
      )
    )
  }


  comp_meta <- obj[[]]


  # Remove cells without replicate/condition/cell-type identity.
  # 不允许缺失关键metadata的cell进入composition。

  valid_comp <- (
    !is.na(
      comp_meta[[BIOLOGICAL_REPLICATE_COL]]
    ) &
      !is.na(
        comp_meta[[CONDITION_COL]]
      ) &
      !is.na(
        comp_meta$cell_type
      ) &
      trimws(
        as.character(
          comp_meta[[BIOLOGICAL_REPLICATE_COL]]
        )
      ) != "" &
      trimws(
        as.character(
          comp_meta[[CONDITION_COL]]
        )
      ) != "" &
      trimws(
        as.character(
          comp_meta$cell_type
        )
      ) != ""
  )


  comp_meta <- comp_meta[
    valid_comp,
    ,
    drop = FALSE
  ]


  if (nrow(comp_meta) == 0L) {

    safe_stop_module(
      MODULE,
      "No valid cells remain for composition analysis."
    )
  }


  # ==========================================================
  # IMPORTANT:
  #
  # group_by includes BOTH biological replicate and condition.
  #
  # 这里同时按replicate + condition分组，
  # 避免理论上相同replicate ID出现在不同condition时
  # 被错误合并。
  # ==========================================================

  composition_cell_type_levels <- sort(
    unique(
      as.character(
        comp_meta$cell_type
      )
    )
  )


  composition <- comp_meta |>
    tibble::rownames_to_column(
      "cell"
    ) |>
    dplyr::count(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]],
      cell_type,
      name = "n_cells"
    )


  names(composition)[
    names(composition) ==
      BIOLOGICAL_REPLICATE_COL
  ] <- ".composition_replicate"


  names(composition)[
    names(composition) ==
      CONDITION_COL
  ] <- ".composition_condition"


  composition <- composition |>
    tidyr::complete(
      tidyr::nesting(
        .composition_replicate,
        .composition_condition
      ),
      cell_type = composition_cell_type_levels,
      fill = list(
        n_cells = 0L
      )
    )


  names(composition)[
    names(composition) ==
      ".composition_replicate"
  ] <- BIOLOGICAL_REPLICATE_COL


  names(composition)[
    names(composition) ==
      ".composition_condition"
  ] <- CONDITION_COL


  composition <- composition |>
    dplyr::group_by(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]]
    ) |>
    dplyr::mutate(
      total_cells = sum(n_cells),
      proportion = n_cells / total_cells
    ) |>
    dplyr::ungroup()


  # Explicit zero rows distinguish an observed absence from a missing sample.
  # This is essential for unbiased replicate-level summaries and group plots.
  # 显式补齐零计数行，可区分“该重复中确实未观察到此细胞类型”和“样本缺失”；
  # 这是进行可靠重复层面汇总和组间作图的重要前提。


  data.table::fwrite(
    composition,
    file.path(
      MDIR,
      "tables",
      "cell_type_composition_by_replicate.csv"
    )
  )


  # Keep the complete table, and combine rare plotting categories as Other.
  # 保留完整组成表；作图时将低于阈值的稀有类别合并为Other。

  composition_plot <- composition |>
    dplyr::mutate(
      cell_type = dplyr::if_else(
        n_cells >= MIN_CELLS_FOR_COMPOSITION,
        as.character(cell_type),
        "Other"
      )
    ) |>
    dplyr::group_by(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]],
      cell_type
    ) |>
    dplyr::summarise(
      n_cells = sum(n_cells),
      .groups = "drop"
    ) |>
    dplyr::filter(
      n_cells > 0L
    ) |>
    dplyr::group_by(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]]
    ) |>
    dplyr::mutate(
      total_cells = sum(n_cells),
      proportion = n_cells / total_cells
    ) |>
    dplyr::ungroup()


  composition_condition_levels <- unique(
    c(
      CONTROL_GROUP,
      CASE_GROUP,
      sort(
        setdiff(
          unique(
            as.character(
              composition_plot[[CONDITION_COL]]
            )
          ),
          c(
            CONTROL_GROUP,
            CASE_GROUP
          )
        )
      )
    )
  )


  composition_replicate_order <- composition_plot |>
    dplyr::distinct(
      .data[[BIOLOGICAL_REPLICATE_COL]],
      .data[[CONDITION_COL]]
    ) |>
    dplyr::mutate(
      condition_order = match(
        as.character(
          .data[[CONDITION_COL]]
        ),
        composition_condition_levels
      )
    ) |>
    dplyr::arrange(
      condition_order,
      .data[[BIOLOGICAL_REPLICATE_COL]]
    ) |>
    dplyr::pull(
      .data[[BIOLOGICAL_REPLICATE_COL]]
    ) |>
    as.character()


  composition_plot_cell_type_levels <- composition_plot |>
    dplyr::group_by(cell_type) |>
    dplyr::summarise(
      total_cells = sum(n_cells),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      dplyr::desc(total_cells),
      cell_type
    ) |>
    dplyr::pull(cell_type) |>
    as.character()


  if (
    "Other" %in%
      composition_plot_cell_type_levels
  ) {
    composition_plot_cell_type_levels <- c(
      setdiff(
        composition_plot_cell_type_levels,
        "Other"
      ),
      "Other"
    )
  }


  composition_plot[[BIOLOGICAL_REPLICATE_COL]] <- factor(
    as.character(
      composition_plot[[BIOLOGICAL_REPLICATE_COL]]
    ),
    levels = rev(
      composition_replicate_order
    )
  )


  composition_plot[[CONDITION_COL]] <- factor(
    as.character(
      composition_plot[[CONDITION_COL]]
    ),
    levels = composition_condition_levels
  )


  composition_plot$cell_type <- factor(
    as.character(
      composition_plot$cell_type
    ),
    levels = composition_plot_cell_type_levels
  )


  composition_cell_type_colours <- stats::setNames(
    grDevices::hcl.colors(
      length(
        composition_plot_cell_type_levels
      ),
      palette = "Dynamic"
    ),
    composition_plot_cell_type_levels
  )


  if (
    "Other" %in%
      names(composition_cell_type_colours)
  ) {
    composition_cell_type_colours[
      "Other"
    ] <- "#BDBDBD"
  }


  # Horizontal 100% stacked bars follow the reference figure while preserving
  # replicate-level proportions and separating conditions into explicit facets.
  # 横向100%堆叠柱图参考文章的展示方式，同时保留每个biological replicate的
  # 独立比例，并使用明确的condition分面避免样本分组含义不清。

  p_comp <- ggplot(
    composition_plot,
    aes(
      x = proportion,
      y = .data[[BIOLOGICAL_REPLICATE_COL]],
      fill = cell_type
    )
  ) +
    geom_col(
      width = 0.82,
      colour = "white",
      linewidth = 0.25
    ) +
    facet_grid(
      stats::as.formula(
        paste0(
          CONDITION_COL,
          "~."
        )
      ),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    scale_x_continuous(
      labels = scales::percent_format(
        accuracy = 1
      ),
      expand = ggplot2::expansion(
        mult = c(
          0,
          0.015
        )
      )
    ) +
    scale_fill_manual(
      values = composition_cell_type_colours,
      drop = FALSE
    ) +
    theme_classic() +
    theme(
      axis.ticks.y = element_blank(),
      strip.placement = "outside",
      strip.background = element_blank(),
      strip.text.y.left = element_text(
        angle = 0,
        face = "bold"
      ),
      plot.title = element_text(
        face = "bold"
      ),
      legend.position = "right",
      plot.margin = ggplot2::margin(
        t = 10,
        r = 15,
        b = 10,
        l = 20
      )
    ) +
    labs(
      title = "Cell-type composition by biological replicate",
      subtitle = "Horizontal bars sum to 100% within each biological replicate",
      x = "Cell proportion",
      y = NULL,
      fill = "Cell type"
    ) +
    guides(
      fill = guide_legend(
        reverse = TRUE
      )
    )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "cell_type_composition.png"
    ),
    p_comp,
    width = 12,
    height = max(
      6,
      0.55 * length(
        composition_replicate_order
      ) + 2.5
    ),
    dpi = 300
  )


  # ==========================================================
  # EXPLORATORY BATCH COMPOSITION COMPARISON
  # 探索性批量细胞组成比较
  # ==========================================================
  #
  # The reference figure performs one group comparison for each cell type.
  # Here, every point is an independent biological replicate, not an individual
  # cell. A Wilcoxon rank-sum test is used because replicate-level proportions
  # are bounded and may not satisfy a Gaussian assumption. P values are adjusted
  # across cell types with the configured method.
  #
  # 参考文章对每种细胞类型分别进行组间比较。本流程中每个点代表一个独立
  # biological replicate，而不是单个cell。由于重复层面的比例数据具有0-1边界，
  # 且未必满足正态分布，这里使用Wilcoxon秩和检验，并按配置方法在所有细胞类型
  # 之间进行多重检验校正。
  #
  # This remains an exploratory screen. Cell-type proportions are compositional
  # and mutually dependent; confirmatory differential-abundance claims require a
  # dedicated model and adequate biological replication.
  #
  # 本结果仍属于探索性筛查。不同细胞类型的比例具有组成约束并相互依赖；若要形成
  # 验证性的差异丰度结论，仍需要专门的组成模型和充分的生物学重复。
  # ==========================================================

  if (RUN_EXPLORATORY_COMPOSITION_TEST) {

    composition_test_data <- composition |>
      dplyr::filter(
        as.character(
          .data[[CONDITION_COL]]
        ) %in% c(
          CONTROL_GROUP,
          CASE_GROUP
        )
      )


    composition_test_results <- lapply(
      composition_cell_type_levels,
      function(target_cell_type) {

        target_data <- composition_test_data |>
          dplyr::filter(
            as.character(cell_type) ==
              target_cell_type
          )


        control_values <- target_data |>
          dplyr::filter(
            as.character(
              .data[[CONDITION_COL]]
            ) == CONTROL_GROUP
          ) |>
          dplyr::pull(proportion)


        case_values <- target_data |>
          dplyr::filter(
            as.character(
              .data[[CONDITION_COL]]
            ) == CASE_GROUP
          ) |>
          dplyr::pull(proportion)


        n_control <- length(control_values)
        n_case <- length(case_values)


        test_status <- if (
          n_control <
            COMPOSITION_TEST_MIN_REPLICATES_PER_GROUP
        ) {
          paste0(
            "insufficient_",
            CONTROL_GROUP,
            "_replicates"
          )
        } else if (
          n_case <
            COMPOSITION_TEST_MIN_REPLICATES_PER_GROUP
        ) {
          paste0(
            "insufficient_",
            CASE_GROUP,
            "_replicates"
          )
        } else {
          "PASS"
        }


        p_value <- NA_real_


        if (identical(test_status, "PASS")) {

          combined_values <- c(
            control_values,
            case_values
          )


          if (
            length(
              unique(combined_values)
            ) == 1L
          ) {
            p_value <- 1
          } else {
            p_value <- tryCatch(
              suppressWarnings(
                stats::wilcox.test(
                  case_values,
                  control_values,
                  exact = FALSE
                )$p.value
              ),
              error = function(e) {
                NA_real_
              }
            )
          }
        }


        data.frame(
          cell_type = target_cell_type,
          comparison = paste0(
            CASE_GROUP,
            "_vs_",
            CONTROL_GROUP
          ),
          test = "Wilcoxon rank-sum",
          n_control = n_control,
          n_case = n_case,
          mean_control = if (
            n_control > 0L
          ) {
            mean(control_values)
          } else {
            NA_real_
          },
          mean_case = if (
            n_case > 0L
          ) {
            mean(case_values)
          } else {
            NA_real_
          },
          median_control = if (
            n_control > 0L
          ) {
            stats::median(control_values)
          } else {
            NA_real_
          },
          median_case = if (
            n_case > 0L
          ) {
            stats::median(case_values)
          } else {
            NA_real_
          },
          mean_difference_case_minus_control = if (
            n_control > 0L &&
              n_case > 0L
          ) {
            mean(case_values) -
              mean(control_values)
          } else {
            NA_real_
          },
          p_value = p_value,
          status = test_status,
          stringsAsFactors = FALSE
        )
      }
    )


    composition_statistics <- dplyr::bind_rows(
      composition_test_results
    ) |>
      dplyr::mutate(
        p_adjusted = stats::p.adjust(
          p_value,
          method = COMPOSITION_P_ADJUST_METHOD
        ),
        significance = dplyr::case_when(
          is.na(p_adjusted) ~ "not tested",
          p_adjusted < 0.001 ~ "***",
          p_adjusted < 0.01 ~ "**",
          p_adjusted < 0.05 ~ "*",
          TRUE ~ "ns"
        )
      ) |>
      dplyr::arrange(cell_type)


    data.table::fwrite(
      composition_statistics,
      file.path(
        MDIR,
        "tables",
        paste0(
          "cell_type_composition_group_",
          "comparison_exploratory.csv"
        )
      )
    )


    composition_test_data$cell_type <- factor(
      as.character(
        composition_test_data$cell_type
      ),
      levels = composition_cell_type_levels
    )


    composition_test_data[[CONDITION_COL]] <- factor(
      as.character(
        composition_test_data[[CONDITION_COL]]
      ),
      levels = c(
        CONTROL_GROUP,
        CASE_GROUP
      )
    )


    composition_annotation <- composition_test_data |>
      dplyr::group_by(cell_type) |>
      dplyr::summarise(
        maximum_proportion = max(proportion),
        .groups = "drop"
      ) |>
      dplyr::left_join(
        composition_statistics,
        by = "cell_type"
      ) |>
      dplyr::mutate(
        panel_scale = pmax(
          maximum_proportion,
          0.02
        ),
        bracket_y = maximum_proportion +
          0.10 * panel_scale,
        bracket_tick_y = maximum_proportion +
          0.055 * panel_scale,
        label_y = maximum_proportion +
          0.18 * panel_scale,
        x_start = CONTROL_GROUP,
        x_end = CASE_GROUP,
        p_label = dplyr::if_else(
          is.finite(p_adjusted),
          paste0(
            "FDR = ",
            format.pval(
              p_adjusted,
              digits = 2,
              eps = 0.001
            ),
            " ",
            significance
          ),
          "not tested"
        )
      )


    composition_condition_colours <- stats::setNames(
      c(
        "#4575B4",
        "#D73027"
      ),
      c(
        CONTROL_GROUP,
        CASE_GROUP
      )
    )


    p_composition_statistics <- ggplot2::ggplot(
      composition_test_data,
      ggplot2::aes(
        x = .data[[CONDITION_COL]],
        y = proportion,
        fill = .data[[CONDITION_COL]]
      )
    ) +
      ggplot2::geom_jitter(
        shape = 21,
        width = 0.10,
        height = 0,
        size = 2.8,
        colour = "#333333",
        stroke = 0.35,
        alpha = 0.90
      ) +
      ggplot2::stat_summary(
        fun.data = ggplot2::mean_se,
        geom = "errorbar",
        width = 0.18,
        colour = "#666666",
        linewidth = 0.55
      ) +
      ggplot2::stat_summary(
        fun = mean,
        geom = "point",
        shape = 23,
        size = 3.0,
        fill = "white",
        colour = "#333333",
        stroke = 0.55
      ) +
      ggplot2::geom_segment(
        data = composition_annotation,
        mapping = ggplot2::aes(
          x = x_start,
          xend = x_end,
          y = bracket_y,
          yend = bracket_y
        ),
        inherit.aes = FALSE,
        colour = "#555555",
        linewidth = 0.45
      ) +
      ggplot2::geom_segment(
        data = composition_annotation,
        mapping = ggplot2::aes(
          x = x_start,
          xend = x_start,
          y = bracket_y,
          yend = bracket_tick_y
        ),
        inherit.aes = FALSE,
        colour = "#555555",
        linewidth = 0.45
      ) +
      ggplot2::geom_segment(
        data = composition_annotation,
        mapping = ggplot2::aes(
          x = x_end,
          xend = x_end,
          y = bracket_y,
          yend = bracket_tick_y
        ),
        inherit.aes = FALSE,
        colour = "#555555",
        linewidth = 0.45
      ) +
      ggplot2::geom_text(
        data = composition_annotation,
        mapping = ggplot2::aes(
          x = x_start,
          y = label_y,
          label = p_label
        ),
        inherit.aes = FALSE,
        position = ggplot2::position_nudge(
          x = 0.5
        ),
        size = 3.0,
        colour = "#333333"
      ) +
      ggplot2::scale_fill_manual(
        values = composition_condition_colours,
        drop = FALSE
      ) +
      ggplot2::scale_y_continuous(
        labels = scales::percent_format(
          accuracy = 0.1
        ),
        expand = ggplot2::expansion(
          mult = c(
            0.05,
            0.24
          )
        )
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(cell_type),
        ncol = 4,
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste0(
          "Exploratory cell-type proportion comparison: ",
          CASE_GROUP,
          " vs ",
          CONTROL_GROUP
        ),
        subtitle = paste0(
          "Points are biological replicates; diamond = mean; ",
          "error bar = mean +/- SEM; Wilcoxon FDR (",
          COMPOSITION_P_ADJUST_METHOD,
          ")"
        ),
        x = NULL,
        y = "Cell proportion",
        fill = "Condition"
      ) +
      ggplot2::theme_classic(
        base_size = 10
      ) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        axis.text.x = ggplot2::element_text(
          angle = 25,
          hjust = 1
        ),
        legend.position = "top",
        plot.margin = ggplot2::margin(
          t = 10,
          r = 15,
          b = 10,
          l = 24
        )
      )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        paste0(
          "cell_type_composition_group_",
          "comparison_exploratory.png"
        )
      ),
      p_composition_statistics,
      width = 16,
      height = max(
        10,
        3.2 * ceiling(
          length(
            composition_cell_type_levels
          ) / 4
        )
      ),
      dpi = 300
    )
  }


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M13. PSEUDOBULK DIFFERENTIAL EXPRESSION / Pseudobulk差异表达
# ============================================================
#
# 【Module scope / 模块范围】
# For each cell type, raw counts are aggregated within independent biological
# replicates and analysed with edgeR using the declared condition contrast and
# optional sample-level covariates. Each cell type is modelled separately.
# 对每种细胞类型，本模块在独立生物学重复内聚合原始计数，并使用edgeR按照已声明的
# condition比较和可选样本级协变量进行分析；每种细胞类型单独建模。
#
# 【Quality gate / 质量门槛】
# A model is fitted only when retained pseudobulks meet cell and replicate minima,
# metadata are unique within replicate, the design matrix is full rank and the
# condition effect is estimable. Skipped cell types and reasons must be reported,
# not silently omitted.
# 仅当保留的pseudobulk满足最少细胞数和重复数、重复内metadata唯一、设计矩阵满秩且
# condition效应可估计时才拟合模型。被跳过的细胞类型及原因必须报告，不能静默省略。
#
# 【Interpretation boundary / 解读边界】
# Effect size, uncertainty, false-discovery rate, expression level and cross-sample
# consistency should be interpreted together. Statistical significance alone does
# not establish biological importance or causality.
# 应综合解释效应量、不确定性、错误发现率、表达水平和跨样本一致性。仅有统计显著性
# 不能证明生物学重要性或因果关系。
#
# 【核心统计原则 / Core statistical principle】
#
# condition-level differential expression的独立统计单位
# 是biological replicate，而不是cell。
#
# The independent statistical unit is the biological replicate,
# NOT individual cells.
#
#
# 每个cell type分别进行：
#
# raw cell-level counts
#       ↓
# biological replicate × cell type aggregation
#       ↓
# pseudobulk count matrix
#       ↓
# edgeR
#       ↓
# Case vs Control
#
#
# 【为什么使用raw counts？】
#
# edgeR需要raw integer-like counts。
#
# 因此这里始终使用RNA counts，
# 不使用M09 logcounts，
# 更不会使用M10B integrated expression。
#
#
# 【运行前检查】
#
# 1. 每个pseudobulk至少MIN_CELLS_PER_PSEUDOBULK cells
# 2. 每组至少MIN_REPLICATES_PER_GROUP replicates
# 3. replicate的condition/covariates必须唯一
# 4. design matrix必须full rank
# 5. condition不能与covariates完全confounded
#
#
# 【输出】
#
# 每个cell type：
# *_edgeR.csv
# *_significant_DEGs.csv
#
# 总结果：
# all_cell_types_edgeR.csv
# all_FDR_significant_DEGs.csv
# all_significant_DEGs.csv
# significant_DEG_counts_by_cell_type.csv
# pseudobulk_logCPM_by_cell_type.rds
# pseudobulk_sample_metadata_by_cell_type.rds
#
# all_FDR_significant_DEGs.csv仅使用FDR阈值，完整保留统计学显著基因；
# all_significant_DEGs.csv进一步要求绝对log2 fold change达到配置阈值，
# 用于更严格的生物学汇报和绘图。
# all_FDR_significant_DEGs.csv retains every statistically significant gene
# using the FDR threshold alone. all_significant_DEGs.csv additionally requires
# the configured absolute log2-fold-change threshold for stricter reporting.
#
# QC：
# pseudobulk_DE_status.csv
#
# Figures / 图片：
# pseudobulk_DE_volcano_facets.png
# pseudobulk_DE_across_cell_types.png
#
# ============================================================


if (RUN_M13_PSEUDOBULK_DE) {

  MODULE <- "M13_pseudobulk_DE"
  MDIR <- module_dir(MODULE)


  obj <- load_checkpoint(
    "M11_annotated"
  )

  validate_complete_annotation(
    obj,
    MODULE
  )


  DefaultAssay(
    obj
  ) <- "RNA"


  # ==========================================================
  # 1. JOIN RAW RNA COUNTS LAYERS IF NECESSARY
  # ==========================================================

  if (
    inherits(
      obj[["RNA"]],
      "Assay5"
    )
  ) {

    obj <- join_assay_layers_if_needed(
      obj,
      assay = "RNA",
      verbose = TRUE
    )
  }


  # ==========================================================
  # 2. EXTRACT RAW COUNTS
  # ==========================================================

  raw_counts <- SeuratObject::LayerData(
    obj,
    assay = "RNA",
    layer = "counts"
  )

  count_values <- if (
    inherits(raw_counts, "sparseMatrix")
  ) {
    methods::slot(
      raw_counts,
      "x"
    )
  } else {
    as.numeric(raw_counts)
  }

  invalid_counts <- (
    anyNA(count_values) ||
      any(!is.finite(count_values)) ||
      any(count_values < 0) ||
      any(
        abs(
          count_values - round(count_values)
        ) > 1e-8
      )
  )

  if (invalid_counts) {
    safe_stop_module(
      MODULE,
      paste0(
        "Pseudobulk edgeR requires non-negative integer raw counts. ",
        "Pseudobulk edgeR必须使用非负整数raw counts。"
      )
    )
  }


  if (
    ncol(raw_counts) !=
    ncol(obj)
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Raw-count columns do not match Seurat cells. ",
        "raw counts列数与Seurat cell数量不一致。"
      )
    )
  }


  # ==========================================================
  # 3. ALIGN CELL METADATA TO RAW COUNTS
  # ==========================================================
  #
  # This is essential.
  #
  # pseudobulk aggregation assumes that each count-matrix
  # column corresponds exactly to the metadata row used
  # for that cell.
  #
  # pseudobulk最重要的前提之一：
  # count matrix列顺序必须与metadata cell顺序完全一致。
  # ==========================================================

  cell_meta <- obj[[]]


  if (
    !all(
      colnames(raw_counts) %in%
      rownames(cell_meta)
    )
  ) {

    safe_stop_module(
      MODULE,
      "Some raw-count cells are missing from metadata."
    )
  }


  cell_meta <- cell_meta[
    colnames(raw_counts),
    ,
    drop = FALSE
  ]


  if (
    !identical(
      rownames(cell_meta),
      colnames(raw_counts)
    )
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Cell metadata and raw counts could not be aligned. ",
        "metadata与raw counts无法严格对齐。"
      )
    )
  }


  # ==========================================================
  # 4. REQUIRED METADATA
  # ==========================================================

  required_pb_meta <- unique(
    c(
      PSEUDOBULK_REPLICATE_COL,
      CONDITION_COL,
      "cell_type",
      DE_COVARIATES
    )
  )


  missing_pb_meta <- setdiff(
    required_pb_meta,
    colnames(cell_meta)
  )


  if (length(missing_pb_meta) > 0L) {

    safe_stop_module(
      MODULE,
      paste0(
        "Missing pseudobulk metadata columns: ",
        paste(
          missing_pb_meta,
          collapse = ", "
        )
      )
    )
  }

  replicate_condition_check <- cell_meta |>
    dplyr::distinct(
      .data[[PSEUDOBULK_REPLICATE_COL]],
      .data[[CONDITION_COL]]
    ) |>
    dplyr::count(
      .data[[PSEUDOBULK_REPLICATE_COL]],
      name = "n_conditions"
    )

  if (any(replicate_condition_check$n_conditions != 1L)) {
    safe_stop_module(
      MODULE,
      paste0(
        "A biological replicate maps to multiple conditions / ",
        "同一biological replicate对应多个condition。"
      )
    )
  }


  # Remove missing core metadata.
  # 缺失replicate/condition/cell_type的cell不能进入DE。

  core_valid <- (
    !is.na(
      cell_meta[[PSEUDOBULK_REPLICATE_COL]]
    ) &
      !is.na(
        cell_meta[[CONDITION_COL]]
      ) &
      !is.na(
        cell_meta$cell_type
      ) &
      trimws(
        as.character(
          cell_meta[[PSEUDOBULK_REPLICATE_COL]]
        )
      ) != "" &
      trimws(
        as.character(
          cell_meta[[CONDITION_COL]]
        )
      ) != "" &
      trimws(
        as.character(
          cell_meta$cell_type
        )
      ) != ""
  )


  if (!all(core_valid)) {

    message(
      sum(!core_valid),
      " cells excluded from pseudobulk because of missing core metadata."
    )
  }


  cell_meta <- cell_meta[
    core_valid,
    ,
    drop = FALSE
  ]


  raw_counts <- raw_counts[
    ,
    rownames(cell_meta),
    drop = FALSE
  ]

  group_valid <- as.character(
    cell_meta[[CONDITION_COL]]
  ) %in% c(
    CONTROL_GROUP,
    CASE_GROUP
  )

  if (!all(group_valid)) {
    message(
      sum(!group_valid),
      " cells excluded because they are outside the requested contrast / ",
      "个cell因不属于指定对比组而排除。"
    )
  }

  cell_meta <- cell_meta[
    group_valid,
    ,
    drop = FALSE
  ]

  raw_counts <- raw_counts[
    ,
    rownames(cell_meta),
    drop = FALSE
  ]


  # ==========================================================
  # 5. CONDITION VALIDATION
  # ==========================================================

  observed_conditions <- unique(
    as.character(
      cell_meta[[CONDITION_COL]]
    )
  )


  if (
    !CONTROL_GROUP %in%
    observed_conditions
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "CONTROL_GROUP not found: ",
        CONTROL_GROUP
      )
    )
  }


  if (
    !CASE_GROUP %in%
    observed_conditions
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "CASE_GROUP not found: ",
        CASE_GROUP
      )
    )
  }


  # ==========================================================
  # 6. LOOP THROUGH CELL TYPES
  # ==========================================================

  all_celltypes <- sort(
    unique(
      as.character(
        cell_meta$cell_type
      )
    )
  )


  de_list <- list()
  status_list <- list()
  pseudobulk_logcpm_list <- list()
  pseudobulk_sample_metadata_list <- list()

  safe_celltype_names <- make.unique(
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      all_celltypes
    ),
    sep = "_"
  )

  names(safe_celltype_names) <- all_celltypes


  for (
    target_cell_type in
    all_celltypes
  ) {


    safe_name <- safe_celltype_names[[target_cell_type]]


    status <- "PASS"
    reason <- ""


    # --------------------------------------------------------
    # 6.1 Cells belonging to this cell type
    # --------------------------------------------------------

    ct_cells <- rownames(cell_meta)[
      as.character(
        cell_meta$cell_type
      ) ==
        target_cell_type
    ]


    ct_meta <- cell_meta[
      ct_cells,
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # 6.2 Count cells per replicate × condition
    # --------------------------------------------------------

    rep_cells <- ct_meta |>
      tibble::rownames_to_column(
        "cell"
      ) |>
      dplyr::count(
        .data[[PSEUDOBULK_REPLICATE_COL]],
        .data[[CONDITION_COL]],
        name = "n_cells"
      )


    eligible_reps <- rep_cells[
      rep_cells$n_cells >=
        MIN_CELLS_PER_PSEUDOBULK,
      ,
      drop = FALSE
    ]


    # --------------------------------------------------------
    # 6.3 Count eligible biological replicates
    # --------------------------------------------------------

    control_rep_ids <- unique(
      as.character(
        eligible_reps[
          eligible_reps[[CONDITION_COL]] ==
            CONTROL_GROUP,
          PSEUDOBULK_REPLICATE_COL,
          drop = TRUE
        ]
      )
    )


    case_rep_ids <- unique(
      as.character(
        eligible_reps[
          eligible_reps[[CONDITION_COL]] ==
            CASE_GROUP,
          PSEUDOBULK_REPLICATE_COL,
          drop = TRUE
        ]
      )
    )


    n_control <- length(
      control_rep_ids
    )


    n_case <- length(
      case_rep_ids
    )


    if (
      n_control <
      MIN_REPLICATES_PER_GROUP
    ) {

      status <- "SKIPPED"

      reason <- paste0(
        "insufficient_control_replicates: ",
        n_control,
        " < ",
        MIN_REPLICATES_PER_GROUP
      )
    }


    if (
      status == "PASS" &&
      n_case <
      MIN_REPLICATES_PER_GROUP
    ) {

      status <- "SKIPPED"

      reason <- paste0(
        "insufficient_case_replicates: ",
        n_case,
        " < ",
        MIN_REPLICATES_PER_GROUP
      )
    }


    # --------------------------------------------------------
    # 6.4 Keep only eligible replicates
    # --------------------------------------------------------

    eligible_rep_ids <- unique(as.character(
        eligible_reps[[PSEUDOBULK_REPLICATE_COL]]
      )
    )


    keep_cells <- rownames(ct_meta)[as.character(
        ct_meta[[PSEUDOBULK_REPLICATE_COL]]
      ) %in%
        eligible_rep_ids
    ]


    if (
      status == "PASS" &&
      length(keep_cells) == 0L
    ) {

      status <- "SKIPPED"

      reason <-
        "no_cells_after_pseudobulk_cell_threshold"
    }


    # ========================================================
    # 7. BUILD REPLICATE METADATA
    # ========================================================

    if (
      status == "PASS"
    ) {

      meta_use <- cell_meta[
        keep_cells,
        ,
        drop = FALSE
      ]


      rep_meta_cols <- unique(
        c(
          PSEUDOBULK_REPLICATE_COL,
          CONDITION_COL,
          DE_COVARIATES
        )
      )


      rep_meta <- meta_use |>
        dplyr::select(
          dplyr::all_of(
            rep_meta_cols
          )
        ) |>
        dplyr::distinct()


      # ------------------------------------------------------
      # Every replicate must correspond to exactly one
      # condition/covariate combination.
      #
      # 每个replicate必须只有唯一condition/covariate组合。
      # ------------------------------------------------------

      rep_id_table <- table(rep_meta[[PSEUDOBULK_REPLICATE_COL]])


      if (
        any(rep_id_table != 1L)
      ) {

        status <- "SKIPPED"

        reason <-
          "replicate_has_inconsistent_condition_or_covariates"
      }
    }


    # ========================================================
    # 8. AGGREGATE RAW COUNTS
    # ========================================================

    if (
      status == "PASS"
    ) {

      # Keep exact cell order used in count matrix.
      # 保证metadata顺序与count matrix完全一致。

      meta_use <- meta_use[
        keep_cells,
        ,
        drop = FALSE
      ]


      rep_factor <- factor(as.character(
          meta_use[[PSEUDOBULK_REPLICATE_COL]]
          ),
        levels = unique(as.character(meta_use[[PSEUDOBULK_REPLICATE_COL]]
          )
        )
      )


      # One-hot cell → replicate design matrix.
      # 每个cell映射到所属biological replicate。

      aggregation_matrix <-
        Matrix::sparse.model.matrix(
          ~ 0 + rep_factor
        )


      pb_counts <- raw_counts[
        ,
        keep_cells,
        drop = FALSE
      ] %*%
        aggregation_matrix


      colnames(
        pb_counts
      ) <- levels(
        rep_factor
      )


      # Reorder replicate metadata to pseudobulk matrix.
      # 将replicate metadata顺序调整为与pb_counts一致。

      rep_meta <- rep_meta[
        match(
          colnames(pb_counts),
          as.character(
            rep_meta[[PSEUDOBULK_REPLICATE_COL]]
            )
          ),
        ,
        drop = FALSE
      ]


      if (
        anyNA(rep_meta[[PSEUDOBULK_REPLICATE_COL]]
        )
      ) {

        status <- "SKIPPED"

        reason <-
          "failed_to_align_pseudobulk_metadata"
      }
    }


    # ========================================================
    # 9. CONDITION AND COVARIATE VALIDATION
    # ========================================================

    if (
      status == "PASS"
    ) {

      rep_meta[[CONDITION_COL]] <- factor(
        as.character(
          rep_meta[[CONDITION_COL]]),
        levels = c(
          CONTROL_GROUP,
          CASE_GROUP
        )
      )


      if (anyNA(
          rep_meta[[CONDITION_COL]]
        )
      ) {

        status <- "SKIPPED"

        reason <-
          "unexpected_condition_value"
      }
    }


    if (
      status == "PASS" &&
      length(DE_COVARIATES) > 0L
    ) {

      if (
        any(
          !stats::complete.cases(
            rep_meta[
              ,
              DE_COVARIATES,
              drop = FALSE
            ]
          )
        )
      ) {

        status <- "SKIPPED"

        reason <-
          "missing_DE_covariate_values"
      }
    }


    # ========================================================
    # 10. DESIGN MATRIX
    # ========================================================

    if (
      status == "PASS"
    ) {

      # Use reformulate() rather than manually constructing
      # formula text.
      #
      # 使用reformulate避免metadata列名中存在特殊字符时
      # 手工paste公式失败。

      design_terms <- c(
        CONDITION_COL,
        DE_COVARIATES
      )


      design_formula <- stats::reformulate(
        design_terms,
        response = NULL,
        intercept = FALSE
      )


      design <- stats::model.matrix(
        design_formula,
        data = rep_meta
      )


      # Full-rank design is essential.
      # design不满秩通常意味着condition与batch/covariate完全混淆。

      if (
        qr(design)$rank <
        ncol(design)
      ) {

        status <- "SKIPPED"

        reason <-
          "design_matrix_not_full_rank_or_condition_confounding"
      }

      if (
        status == "PASS" &&
        nrow(design) <= ncol(design)
      ) {
        status <- "SKIPPED"
        reason <- "no_residual_degrees_of_freedom"
      }
    }


    # ========================================================
    # 11. edgeR FILTERING
    # ========================================================

    if (
      status == "PASS"
    ) {

      y <- edgeR::DGEList(
        counts = pb_counts,
        samples = rep_meta
      )


      keep_gene <- edgeR::filterByExpr(
        y,
        design = design
      )


      if (
        sum(keep_gene) < 2L
      ) {

        status <- "SKIPPED"

        reason <-
          "too_few_expressed_genes_after_filterByExpr"
      }
    }


    # ========================================================
    # 12. edgeR MODEL
    # ========================================================

    if (
      status == "PASS"
    ) {

      y <- y[
        keep_gene,
        ,
        keep.lib.sizes = FALSE
      ]


      # edgeR renamed calcNormFactors() to normLibSizes(). Use the current name
      # when available and retain a compatibility fallback for older releases.
      # The calculation remains TMM library-size normalization in both cases.
      # edgeR已将calcNormFactors()更名为normLibSizes()。优先使用当前函数名，
      # 同时为旧版edgeR保留兼容分支；两者执行的均为TMM文库量标准化。
      if (
        "normLibSizes" %in%
          getNamespaceExports("edgeR")
      ) {
        y <- edgeR::normLibSizes(
          y
        )
      } else {
        y <- edgeR::calcNormFactors(
          y
        )
      }


      y <- edgeR::estimateDisp(
        y,
        design
      )


      fit <- edgeR::glmQLFit(
        y,
        design,
        robust = TRUE
      )


      # ------------------------------------------------------
      # Identify design columns using model.matrix metadata
      # rather than fragile grep whenever possible.
      # ------------------------------------------------------

      control_name <- paste0(
        CONDITION_COL,
        CONTROL_GROUP
      )


      case_name <- paste0(
        CONDITION_COL,
        CASE_GROUP
      )


      if (
        !control_name %in%
        colnames(design) ||
        !case_name %in%
        colnames(design)
      ) {

        status <- "SKIPPED"

        reason <-
          "cannot_identify_control_case_design_columns"
      }
    }


    # ========================================================
    # 13. CASE VS CONTROL CONTRAST
    # ========================================================

    if (
      status == "PASS"
    ) {

      contrast <- rep(
        0,
        ncol(design)
      )


      names(
        contrast
      ) <- colnames(
        design
      )


      contrast[
        case_name
      ] <- 1


      contrast[
        control_name
      ] <- -1


      qlf <- edgeR::glmQLFTest(
        fit,
        contrast = contrast
      )


      result <- edgeR::topTags(
        qlf,
        n = Inf,
        sort.by = "PValue"
      )$table


      result <- result |>
        tibble::rownames_to_column(
          "gene"
        ) |>
        dplyr::mutate(
          cell_type =
            target_cell_type,
          contrast =
            paste0(
              CASE_GROUP,
              "_vs_",
              CONTROL_GROUP
            ),
          n_control =
            n_control,
          n_case =
            n_case
        ) |>
        dplyr::relocate(
          cell_type,
          contrast,
          gene
        )


      data.table::fwrite(
        result,
        file.path(
          MDIR,
          "tables",
          paste0(
            safe_name,
            "_",
            CASE_GROUP,
            "_vs_",
            CONTROL_GROUP,
            "_edgeR.csv"
          )
        )
      )


      # ------------------------------------------------------
      # Save a thresholded DEG table for convenient biological
      # review while retaining the complete edgeR table above.
      # ------------------------------------------------------
      #
      # The statistical test is not repeated here. Genes are
      # classified from the edgeR FDR and logFC columns using
      # the reporting thresholds declared in the configuration.
      #
      # 在保留上方完整edgeR结果的同时，另存一份便于生物学审阅的
      # 显著DEG表。这里不重复进行统计检验，只依据配置区定义的
      # FDR和logFC阈值对edgeR结果进行分类。
      # ------------------------------------------------------

      significant_result <- result |>
        dplyr::filter(
          is.finite(FDR),
          FDR <= PSEUDOBULK_DE_FDR_THRESHOLD,
          abs(logFC) >=
            PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ) |>
        dplyr::mutate(
          direction = ifelse(
            logFC > 0,
            paste0(
              "Up in ",
              CASE_GROUP
            ),
            paste0(
              "Up in ",
              CONTROL_GROUP
            )
          )
        ) |>
        dplyr::arrange(
          FDR,
          dplyr::desc(
            abs(logFC)
          )
        )


      data.table::fwrite(
        significant_result,
        file.path(
          MDIR,
          "tables",
          paste0(
            safe_name,
            "_",
            CASE_GROUP,
            "_vs_",
            CONTROL_GROUP,
            "_significant_DEGs.csv"
          )
        )
      )


      de_list[[target_cell_type]] <- result


      # ------------------------------------------------------
      # Retain the exact filtered, TMM-normalized pseudobulk
      # expression matrix and aligned replicate metadata used
      # by the valid M13 model. M13C reads these objects so GSVA
      # uses the same biological replicates and expression filter
      # instead of rebuilding a subtly different input.
      #
      # 保留当前通过M13质量门槛的、经过表达过滤及TMM
      # 标准化的pseudobulk表达矩阵，并同时保留严格对齐
      # 的replicate metadata。M13C直接读取这些对象，从而确保
      # GSVA与edgeR使用相同生物学重复和表达过滤。
      # ------------------------------------------------------

      pseudobulk_logcpm_list[[target_cell_type]] <- edgeR::cpm(
        y,
        log = TRUE,
        prior.count = GSVA_LOGCPM_PRIOR_COUNT
      )


      pseudobulk_sample_metadata_list[[target_cell_type]] <- rep_meta |>
        dplyr::mutate(
          pseudobulk_id = as.character(
            .data[[PSEUDOBULK_REPLICATE_COL]]
          ),
          cell_type = target_cell_type,
          .before = 1
        )
    }


    # ========================================================
    # 14. RECORD STATUS FOR THIS CELL TYPE
    # ========================================================

    status_list[[target_cell_type]] <- data.frame(

      cell_type =
        target_cell_type,

      status =
        status,

      reason =
        reason,

      n_control_eligible =
        n_control,

      n_case_eligible =
        n_case,

      min_cells_per_pseudobulk =
        MIN_CELLS_PER_PSEUDOBULK,

      min_replicates_per_group =
        MIN_REPLICATES_PER_GROUP,

      stringsAsFactors = FALSE
    )
  }


  # ==========================================================
  # 15. SAVE STATUS TABLE
  # ==========================================================

  de_status <- dplyr::bind_rows(
    status_list
  )


  data.table::fwrite(
    de_status,
    file.path(
      MDIR,
      "tables",
      "pseudobulk_DE_status.csv"
    )
  )


  # ==========================================================
  # 16. COMBINE VALID RESULTS
  # ==========================================================

  if (
    length(de_list) > 0L
  ) {

    saveRDS(
      pseudobulk_logcpm_list,
      file.path(
        MDIR,
        "tables",
        "pseudobulk_logCPM_by_cell_type.rds"
      )
    )


    saveRDS(
      pseudobulk_sample_metadata_list,
      file.path(
        MDIR,
        "tables",
        "pseudobulk_sample_metadata_by_cell_type.rds"
      )
    )

    de_all <- dplyr::bind_rows(
      de_list
    )


    data.table::fwrite(
      de_all,
      file.path(
        MDIR,
        "tables",
        "all_cell_types_edgeR.csv"
      )
    )


    # ========================================================
    # 17. CLASSIFY AND SUMMARIZE SIGNIFICANT DEGs
    # ========================================================
    #
    # edgeR has already tested every retained gene separately
    # within each cell type. This section combines the valid
    # results for reporting; it does not pool cells or refit a
    # model across cell types.
    #
    # edgeR已在每种细胞类型内部对保留基因完成独立检验。本节仅将
    # 通过质量门槛的结果合并用于汇报，不会合并不同细胞类型的cell，
    # 也不会跨细胞类型重新拟合统计模型。
    #
    # A gene is classified as significant only when it passes
    # both the FDR and absolute-logFC reporting thresholds.
    # Complete edgeR tables remain available without this filter.
    #
    # 只有同时满足FDR和绝对logFC汇报阈值的基因才被归为显著DEG；
    # 未经该阈值筛选的完整edgeR结果仍会保留。
    # ========================================================

    up_case_label <- paste0(
      "Up in ",
      CASE_GROUP
    )


    up_control_label <- paste0(
      "Up in ",
      CONTROL_GROUP
    )


    de_significance_levels <- c(
      up_case_label,
      up_control_label,
      "Not significant"
    )


    de_cell_type_levels <- all_celltypes[
      all_celltypes %in%
        unique(de_all$cell_type)
    ]


    de_all_visual <- de_all |>
      dplyr::mutate(
        direction = ifelse(
          logFC >= 0,
          up_case_label,
          up_control_label
        ),
        significance = dplyr::if_else(
          is.finite(FDR) &
            FDR <= PSEUDOBULK_DE_FDR_THRESHOLD &
            abs(logFC) >=
              PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD,
          direction,
          "Not significant"
        ),
        significance = factor(
          significance,
          levels = de_significance_levels
        ),
        cell_type = factor(
          cell_type,
          levels = de_cell_type_levels
        ),
        negative_log10_FDR = -log10(
          pmax(
            FDR,
            .Machine$double.xmin
          )
        )
      )


    # ========================================================
    # 17A. COMPLETE FDR-SIGNIFICANT GENE TABLE
    #      完整FDR显著基因表
    # ========================================================
    #
    # This table retains every gene meeting the configured FDR
    # threshold, regardless of effect size. The additional
    # passes_abs_logFC_threshold column indicates whether each
    # gene also enters the stricter DEG table used for figures
    # and concise biological reporting.
    #
    # 本表完整保留所有达到配置FDR阈值的基因，不再额外删除效应量
    # 较小的基因。passes_abs_logFC_threshold列用于说明该基因是否
    # 同时达到严格DEG表及绘图所使用的绝对log2 fold change阈值。
    # ========================================================

    de_fdr_significant_all <- de_all_visual |>
      dplyr::filter(
        is.finite(FDR),
        FDR <= PSEUDOBULK_DE_FDR_THRESHOLD
      ) |>
      dplyr::mutate(
        cell_type = as.character(
          cell_type
        ),
        direction = as.character(
          direction
        ),
        passes_abs_logFC_threshold =
          abs(logFC) >=
            PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD,
        reporting_class = dplyr::if_else(
          passes_abs_logFC_threshold,
          "FDR and effect-size thresholds",
          "FDR threshold only"
        )
      ) |>
      dplyr::select(
        -significance,
        -negative_log10_FDR
      ) |>
      dplyr::arrange(
        match(
          cell_type,
          de_cell_type_levels
        ),
        FDR,
        dplyr::desc(
          abs(logFC)
        )
      )


    data.table::fwrite(
      de_fdr_significant_all,
      file.path(
        MDIR,
        "tables",
        "all_FDR_significant_DEGs.csv"
      )
    )


    de_significant_all <- de_all_visual |>
      dplyr::filter(
        significance != "Not significant"
      ) |>
      dplyr::arrange(
        cell_type,
        FDR,
        dplyr::desc(
          abs(logFC)
        )
      )


    data.table::fwrite(
      de_significant_all |>
        dplyr::mutate(
          cell_type = as.character(cell_type),
          significance = as.character(significance)
        ),
      file.path(
        MDIR,
        "tables",
        "all_significant_DEGs.csv"
      )
    )


    significant_de_counts <- expand.grid(
      cell_type = de_cell_type_levels,
      direction = c(
        up_case_label,
        up_control_label
      ),
      stringsAsFactors = FALSE
    ) |>
      dplyr::left_join(
        de_significant_all |>
          dplyr::mutate(
            cell_type = as.character(cell_type)
          ) |>
          dplyr::count(
            cell_type,
            direction,
            name = "n_significant_DEGs"
          ),
        by = c(
          "cell_type",
          "direction"
        )
      ) |>
      dplyr::mutate(
        n_significant_DEGs = tidyr::replace_na(
          n_significant_DEGs,
          0L
        )
      ) |>
      dplyr::arrange(
        match(
          cell_type,
          de_cell_type_levels
        ),
        match(
          direction,
          c(
            up_case_label,
            up_control_label
          )
        )
      )


    data.table::fwrite(
      significant_de_counts,
      file.path(
        MDIR,
        "tables",
        "significant_DEG_counts_by_cell_type.csv"
      )
    )


    # ========================================================
    # 18. SELECT LABELS FOR DIFFERENTIAL-EXPRESSION FIGURES
    # ========================================================
    #
    # Only a limited number of genes are labelled in each cell
    # type and direction. Selection uses the smallest FDR first
    # and the largest absolute logFC second. All tested genes
    # remain visible as points and remain present in the tables.
    #
    # 每种细胞类型和变化方向只标注有限数量的基因。标签首先选择FDR
    # 最小的基因，再依据绝对logFC排序。所有接受检验的基因仍以散点
    # 显示，并完整保留在结果表中。
    # ========================================================

    de_label_data <- de_significant_all |>
      dplyr::group_by(
        cell_type,
        direction
      ) |>
      dplyr::arrange(
        FDR,
        dplyr::desc(
          abs(logFC)
        ),
        .by_group = TRUE
      ) |>
      dplyr::slice_head(
        n = PSEUDOBULK_DE_LABELS_PER_DIRECTION
      ) |>
      dplyr::ungroup()


    de_plot_colours <- c(
      stats::setNames(
        "#D73027",
        up_case_label
      ),
      stats::setNames(
        "#4575B4",
        up_control_label
      ),
      `Not significant` = "#BDBDBD"
    )


    # ========================================================
    # 19. MULTI-CELL-TYPE DIFFERENTIAL-EXPRESSION PLOT
    # ========================================================
    #
    # This figure adapts the article's combined multi-group DEG
    # display to replicate-aware edgeR pseudobulk results. The
    # x-axis separates cell types and the y-axis shows CASE_GROUP
    # versus CONTROL_GROUP log2 fold change. Horizontal dashed
    # lines mark the configured effect-size threshold.
    #
    # 本图将参考文章的多分组DEG合并展示方式改造成适用于edgeR
    # pseudobulk结果的版本。横轴区分细胞类型，纵轴表示CASE_GROUP
    # 相对CONTROL_GROUP的log2 fold change；水平虚线表示配置区定义的
    # 效应量阈值。
    # ========================================================

    set.seed(RANDOM_SEED)


    de_all_visual <- de_all_visual |>
      dplyr::mutate(
        cell_type_index = as.integer(cell_type),
        plot_x = cell_type_index +
          stats::runif(
            dplyr::n(),
            min = -0.28,
            max = 0.28
          )
      )


    de_label_data <- de_label_data |>
      dplyr::left_join(
        de_all_visual |>
          dplyr::select(
            cell_type,
            gene,
            plot_x
          ),
        by = c(
          "cell_type",
          "gene"
        )
      )


    de_panel_background <- data.frame(
      xmin = seq_along(de_cell_type_levels) - 0.5,
      xmax = seq_along(de_cell_type_levels) + 0.5,
      ymin = -Inf,
      ymax = Inf,
      shade = seq_along(de_cell_type_levels) %% 2L == 1L
    )


    p_de_across_cell_types <- ggplot2::ggplot() +
      ggplot2::geom_rect(
        data = de_panel_background[
          de_panel_background$shade,
          ,
          drop = FALSE
        ],
        mapping = ggplot2::aes(
          xmin = xmin,
          xmax = xmax,
          ymin = ymin,
          ymax = ymax
        ),
        inherit.aes = FALSE,
        fill = "#BDBDBD",
        alpha = 0.16,
        colour = NA
      ) +
      ggplot2::geom_hline(
        yintercept = 0,
        colour = "#333333",
        linewidth = 0.45
      ) +
      ggplot2::geom_hline(
        yintercept = c(
          -PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD,
          PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ),
        colour = "#666666",
        linewidth = 0.35,
        linetype = "dashed"
      ) +
      ggplot2::geom_point(
        data = de_all_visual |>
          dplyr::filter(
            significance == "Not significant"
          ),
        mapping = ggplot2::aes(
          x = plot_x,
          y = logFC
        ),
        colour = "#BDBDBD",
        alpha = 0.24,
        size = 0.55
      ) +
      ggplot2::geom_point(
        data = de_all_visual |>
          dplyr::filter(
            significance != "Not significant"
          ),
        mapping = ggplot2::aes(
          x = plot_x,
          y = logFC,
          colour = significance
        ),
        alpha = 0.90,
        size = 1.15
      ) +
      ggrepel::geom_text_repel(
        data = de_label_data,
        mapping = ggplot2::aes(
          x = plot_x,
          y = logFC,
          label = gene,
          colour = significance
        ),
        seed = RANDOM_SEED,
        size = 3.0,
        box.padding = 0.30,
        point.padding = 0.15,
        min.segment.length = 0,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = de_plot_colours,
        drop = FALSE
      ) +
      ggplot2::scale_x_continuous(
        breaks = seq_along(de_cell_type_levels),
        labels = de_cell_type_levels,
        expand = ggplot2::expansion(
          add = 0.65
        )
      ) +
      ggplot2::scale_y_continuous(
        expand = ggplot2::expansion(
          mult = c(
            0.10,
            0.18
          )
        )
      ) +
      ggplot2::labs(
        title = paste0(
          "Pseudobulk differential expression across cell types: ",
          CASE_GROUP,
          " vs ",
          CONTROL_GROUP
        ),
        subtitle = paste0(
          "Significant DEG: FDR <= ",
          PSEUDOBULK_DE_FDR_THRESHOLD,
          " and |log2FC| >= ",
          PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD,
          "; positive log2FC = higher in ",
          CASE_GROUP
        ),
        x = "Cell type",
        y = paste0(
          "log2FC (",
          CASE_GROUP,
          " vs ",
          CONTROL_GROUP,
          ")"
        ),
        colour = "DEG classification"
      ) +
      ggplot2::theme_minimal(
        base_size = 11
      ) +
      ggplot2::theme(
        panel.grid.major.x = ggplot2::element_blank(),
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1,
          vjust = 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "top",
        plot.margin = ggplot2::margin(
          t = 10,
          r = 15,
          b = 10,
          l = 28
        )
      )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "pseudobulk_DE_across_cell_types.png"
      ),
      p_de_across_cell_types,
      width = max(
        16,
        1.2 * length(de_cell_type_levels)
      ),
      height = 10,
      dpi = 300
    )


    # ========================================================
    # 20. FACETED PSEUDOBULK VOLCANO PLOT
    # ========================================================
    #
    # A faceted volcano overview complements the combined plot.
    # Each panel remains a separate cell-type comparison even
    # though all valid panels are displayed in one image.
    #
    # 多面板火山图用于补充合并差异散点图。尽管所有通过质量门槛的
    # 细胞类型展示在同一张图片中，每个panel仍对应独立的细胞类型
    # pseudobulk比较。
    # ========================================================

    p_de_volcano_facets <- ggplot2::ggplot(
      de_all_visual,
      ggplot2::aes(
        x = logFC,
        y = negative_log10_FDR,
        colour = significance
      )
    ) +
      ggplot2::geom_hline(
        yintercept = -log10(
          PSEUDOBULK_DE_FDR_THRESHOLD
        ),
        colour = "#666666",
        linewidth = 0.35,
        linetype = "dashed"
      ) +
      ggplot2::geom_vline(
        xintercept = c(
          -PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD,
          PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ),
        colour = "#666666",
        linewidth = 0.35,
        linetype = "dashed"
      ) +
      ggplot2::geom_point(
        alpha = 0.45,
        size = 0.60
      ) +
      ggrepel::geom_text_repel(
        data = de_label_data,
        ggplot2::aes(
          x = logFC,
          y = negative_log10_FDR,
          label = gene,
          colour = significance
        ),
        seed = RANDOM_SEED,
        size = 2.6,
        box.padding = 0.25,
        point.padding = 0.10,
        min.segment.length = 0,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      ggplot2::scale_colour_manual(
        values = de_plot_colours,
        drop = FALSE
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(cell_type),
        ncol = 4,
        scales = "free"
      ) +
      ggplot2::labs(
        title = paste0(
          "Cell-type pseudobulk volcano plots: ",
          CASE_GROUP,
          " vs ",
          CONTROL_GROUP
        ),
        subtitle = paste0(
          "Dashed thresholds: FDR = ",
          PSEUDOBULK_DE_FDR_THRESHOLD,
          ", |log2FC| = ",
          PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ),
        x = "log2 fold change",
        y = expression(
          -log[10](FDR)
        ),
        colour = "DEG classification"
      ) +
      ggplot2::theme_minimal(
        base_size = 10
      ) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "top",
        plot.margin = ggplot2::margin(
          t = 10,
          r = 15,
          b = 10,
          l = 28
        )
      )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "pseudobulk_DE_volcano_facets.png"
      ),
      p_de_volcano_facets,
      width = 16,
      height = max(
        10,
        3.2 * ceiling(
          length(de_cell_type_levels) / 4
        )
      ),
      dpi = 300
    )

  } else {

    message(
      "No cell type passed pseudobulk DE validation / ",
      "没有cell type通过pseudobulk DE统计条件。"
    )
  }


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M13B. FUNCTIONAL INTERPRETATION OF PSEUDOBULK DE
# M13B. Pseudobulk差异基因功能解释
# ============================================================
#
# 【Module scope / 模块范围】
# This module interprets the replicate-aware M13 results. It does not run a new
# differential-expression test and it does not replace the complete edgeR tables.
# The same configured FDR and absolute-log2FC thresholds are used to define the
# input DEG classes, so tables, heatmaps and enrichment reports remain traceable
# to one explicit statistical decision.
# 本模块对具有生物学重复感知的M13结果进行功能解释，不重新执行差异
# 检验，也不取代完整edgeR表。模块使用配置区中相同的FDR和绝对log2FC
# 门槛定义DEG类别，使表格、热图及富集结果都能追溯到同一明确的统计判定。
#
# 【Statistical universe / 统计背景】
# GO and KEGG over-representation analysis uses, separately for every cell type,
# all genes retained and tested by edgeR as the universe. This is preferable to
# the whole genome because it conditions the interpretation on genes detectable
# under the current experiment and M13 filtering process.
# GO与KEGG过度富集会针对每种cell type，将edgeR实际保留并检验的全部基因
# 作为universe。这比使用全基因组更合理，因为富集解释会以当前实验可检测
# 且通过M13过滤的基因集为条件。
#
# 【Interpretation boundary / 解读边界】
# Enrichment identifies annotation terms represented more often than expected;
# it does not show that a pathway is experimentally activated. ggkegg overlays
# pseudobulk log2FC values on database pathway nodes and therefore remains an
# annotated visualization, not a pathway-activity or causal-signalling test.
# 富集分析反映注释术语在DEG中出现的频率高于背景预期，不能单独证明通路
# 已在实验中激活。ggkegg将pseudobulk log2FC叠加到数据库通路节点上，因此属于
# 带注释的可视化，不是通路活性或因果信号检验。
#
# 【Principal output / 主要输出】
#   tables/DEG_class_heatmap_selection.csv
#   tables/DEG_class_logFC_matrix.csv
#   tables/gene_identifier_mapping.csv
#   tables/enrichment_class_audit.csv
#   tables/GO_BP_enrichment_all.csv
#   tables/KEGG_enrichment_all.csv
#   tables/ggkegg_pathway_audit.csv
#   tables/ggkegg_pathway_selection.csv
#   tables/ggkegg_mapped_significant_DEGs.csv
#   figures/DEG_class_logFC_heatmap.png
#   figures/GO_BP_enrichment_dotplot.png
#   figures/KEGG_enrichment_dotplot.png
#   figures/ggkegg/pathway_maps/*.png
#   figures/ggkegg/composite_panels/*.png
#
# 【Method references / 方法依据】
# clusterProfiler: https://doi.org/10.1089/omi.2011.0118
# clusterProfiler manual:
# https://bioconductor.org/packages/release/bioc/manuals/clusterProfiler/man/clusterProfiler.pdf
# ggkegg package and manual:
# https://bioconductor.org/packages/release/bioc/html/ggkegg.html
# ============================================================


if (RUN_M13B_FUNCTIONAL) {

  MODULE <- "M13B_functional_interpretation"
  MDIR <- module_dir(MODULE)


  required_functional_packages <- c(
    "AnnotationDbi",
    "clusterProfiler",
    "ComplexHeatmap",
    "circlize",
    species_cfg$orgdb_package
  )


  if (
    RUN_GGKEGG
  ) {
    required_functional_packages <- c(
      required_functional_packages,
      "ggkegg",
      "ggraph",
      "tidygraph",
      "igraph",
      "patchwork"
    )
  }


  functional_packages_available <- vapply(
    unique(
      required_functional_packages
    ),
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )


  if (
    any(
      !functional_packages_available
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Missing functional-interpretation packages / ",
        "缺少功能解释所需软件包: ",
        paste(
          names(
            functional_packages_available
          )[
            !functional_packages_available
          ],
          collapse = ", "
        ),
        ". Run M00 with INSTALL_PACKAGES = TRUE."
      )
    )
  }


  m13_table_dir <- file.path(
    DIR_MODULES,
    "M13_pseudobulk_DE",
    "tables"
  )


  complete_de_file <- file.path(
    m13_table_dir,
    "all_cell_types_edgeR.csv"
  )


  significant_de_file <- file.path(
    m13_table_dir,
    "all_significant_DEGs.csv"
  )


  if (
    !file.exists(
      complete_de_file
    ) ||
    !file.exists(
      significant_de_file
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "M13 combined DEG tables are missing / M13合并DEG表不存在。\n",
        "Run M13 first: ",
        m13_table_dir
      )
    )
  }


  de_all <- data.table::fread(
    complete_de_file,
    data.table = FALSE
  )


  de_significant <- data.table::fread(
    significant_de_file,
    data.table = FALSE
  )


  required_de_columns <- c(
    "cell_type",
    "gene",
    "logFC",
    "FDR"
  )


  if (
    any(
      !required_de_columns %in%
        colnames(
          de_all
        )
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "M13 complete table lacks required columns / M13完整表缺少必需列: ",
        paste(
          setdiff(
            required_de_columns,
            colnames(de_all)
          ),
          collapse = ", "
        )
      )
    )
  }


  if (
    nrow(
      de_significant
    ) == 0L
  ) {
    data.table::fwrite(
      data.frame(
        status = "SKIPPED",
        reason = paste0(
          "No gene passed FDR <= ",
          PSEUDOBULK_DE_FDR_THRESHOLD,
          " and |log2FC| >= ",
          PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ),
        stringsAsFactors = FALSE
      ),
      file.path(
        MDIR,
        "tables",
        "functional_interpretation_status.csv"
      )
    )

    write_module_status(
      MODULE,
      "OK",
      "No significant DEG class was available; audit table written."
    )

  } else {

    up_case_label <- paste0(
      "Up in ",
      CASE_GROUP
    )


    up_control_label <- paste0(
      "Up in ",
      CONTROL_GROUP
    )


    de_significant <- de_significant |>
      dplyr::mutate(
        cell_type = as.character(
          .data$cell_type
        ),
        gene = as.character(
          .data$gene
        ),
        direction = dplyr::if_else(
          .data$logFC > 0,
          up_case_label,
          up_control_label
        ),
        class = paste(
          .data$cell_type,
          .data$direction,
          sep = " | "
        )
      )


    de_cell_type_levels <- unique(
      as.character(
        de_all$cell_type
      )
    )


    # ========================================================
    # 1. CLASSIFIED CROSS-CELL-TYPE log2FC HEATMAP
    # ========================================================
    # Rows are selected within the cell type and direction in which a DEG was
    # discovered. Columns show the same gene's pseudobulk log2FC in every tested
    # cell type. Thus the figure distinguishes a source DEG class while revealing
    # whether the direction is shared, absent or reversed elsewhere.
    # 行在DEG被发现的cell type和方向内选取；列展示同一基因在所有已检验
    # cell type中的pseudobulk log2FC。因此该图既保留来源DEG分类，也能显示
    # 相同方向是否在其他cell type中共享、缺失或反转。

    heatmap_selection <- de_significant |>
      dplyr::group_by(
        .data$cell_type,
        .data$direction
      ) |>
      dplyr::arrange(
        .data$FDR,
        dplyr::desc(
          abs(
            .data$logFC
          )
        ),
        .by_group = TRUE
      ) |>
      dplyr::slice_head(
        n = as.integer(
          DEG_CLASS_HEATMAP_TOP_PER_CLASS
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::mutate(
        row_id = make.unique(
          paste(
            .data$gene,
            .data$cell_type,
            .data$direction,
            sep = " | "
          )
        )
      )


    data.table::fwrite(
      heatmap_selection,
      file.path(
        MDIR,
        "tables",
        "DEG_class_heatmap_selection.csv"
      )
    )


    heatmap_long <- tidyr::crossing(
      row_id = heatmap_selection$row_id,
      comparison_cell_type = de_cell_type_levels
    ) |>
      dplyr::left_join(
        heatmap_selection |>
          dplyr::select(
            "row_id",
            source_cell_type = "cell_type",
            source_direction = "direction",
            "gene"
          ),
        by = "row_id"
      ) |>
      dplyr::left_join(
        de_all |>
          dplyr::transmute(
            comparison_cell_type = as.character(
              .data$cell_type
            ),
            gene = as.character(
              .data$gene
            ),
            logFC = as.numeric(
              .data$logFC
            ),
            FDR = as.numeric(
              .data$FDR
            )
          ),
        by = c(
          "comparison_cell_type",
          "gene"
        )
      )


    data.table::fwrite(
      heatmap_long,
      file.path(
        MDIR,
        "tables",
        "DEG_class_logFC_matrix.csv"
      )
    )


    heatmap_matrix <- heatmap_long |>
      dplyr::select(
        "row_id",
        "comparison_cell_type",
        "logFC"
      ) |>
      tidyr::pivot_wider(
        names_from = "comparison_cell_type",
        values_from = "logFC"
      ) |>
      tibble::column_to_rownames(
        "row_id"
      ) |>
      as.matrix()


    heatmap_matrix <- heatmap_matrix[
      heatmap_selection$row_id,
      de_cell_type_levels,
      drop = FALSE
    ]


    row_labels <- stats::setNames(
      heatmap_selection$gene,
      heatmap_selection$row_id
    )


    row_class <- factor(
      paste(
        heatmap_selection$cell_type,
        heatmap_selection$direction,
        sep = " | "
      ),
      levels = unique(
        paste(
          heatmap_selection$cell_type,
          heatmap_selection$direction,
          sep = " | "
        )
      )
    )


    finite_logfc <- abs(
      heatmap_matrix[
        is.finite(
          heatmap_matrix
        )
      ]
    )


    heatmap_limit <- if (
      length(
        finite_logfc
      ) == 0L
    ) {
      1
    } else {
      max(
        1,
        min(
          4,
          stats::quantile(
            finite_logfc,
            probs = 0.98,
            na.rm = TRUE,
            names = FALSE
          )
        )
      )
    }


    de_heatmap <- ComplexHeatmap::Heatmap(
      heatmap_matrix,
      name = "log2FC",
      col = circlize::colorRamp2(
        c(
          -heatmap_limit,
          0,
          heatmap_limit
        ),
        c(
          "#2166AC",
          "white",
          "#B2182B"
        )
      ),
      na_col = "#E6E6E6",
      cluster_rows = FALSE,
      cluster_columns = FALSE,
      row_split = row_class,
      row_labels = unname(
        row_labels[
          rownames(
            heatmap_matrix
          )
        ]
      ),
      show_row_names = TRUE,
      show_column_names = TRUE,
      row_names_gp = grid::gpar(
        fontsize = 7.5
      ),
      column_names_gp = grid::gpar(
        fontsize = 9
      ),
      column_names_rot = 45,
      row_title_gp = grid::gpar(
        fontsize = 8,
        fontface = "bold"
      ),
      row_gap = grid::unit(
        1.5,
        "mm"
      ),
      border = TRUE,
      column_title = paste0(
        "Pseudobulk DEG classes: ",
        CASE_GROUP,
        " vs ",
        CONTROL_GROUP
      ),
      heatmap_legend_param = list(
        title = paste0(
          "log2FC\n",
          CASE_GROUP,
          " vs ",
          CONTROL_GROUP
        )
      )
    )


    ragg::agg_png(
      filename = file.path(
        MDIR,
        "figures",
        "DEG_class_logFC_heatmap.png"
      ),
      width = max(
        11,
        0.65 * length(
          de_cell_type_levels
        ) + 5
      ),
      height = max(
        9,
        0.20 * nrow(
          heatmap_matrix
        ) + 3
      ),
      units = "in",
      res = 300,
      background = "transparent"
    )


    tryCatch(
      ComplexHeatmap::draw(
        de_heatmap,
        background = "transparent",
        heatmap_legend_side = "right"
      ),
      finally = {
        if (
          grDevices::dev.cur() > 1L
        ) {
          grDevices::dev.off()
        }
      }
    )


    # ========================================================
    # 2. SYMBOL-TO-ENTREZ MAPPING
    # ========================================================

    orgdb_object <- getExportedValue(
      species_cfg$orgdb_package,
      species_cfg$orgdb_package
    )


    all_symbols <- sort(
      unique(
        as.character(
          de_all$gene
        )
      )
    )


    identifier_mapping <- suppressMessages(
      AnnotationDbi::select(
        orgdb_object,
        keys = all_symbols,
        keytype = "SYMBOL",
        columns = c(
          "SYMBOL",
          "ENTREZID"
        )
      )
    ) |>
      dplyr::filter(
        !is.na(
          .data$ENTREZID
        ),
        trimws(
          as.character(
            .data$ENTREZID
          )
        ) != ""
      ) |>
      dplyr::distinct(
        .data$SYMBOL,
        .data$ENTREZID
      )


    data.table::fwrite(
      identifier_mapping,
      file.path(
        MDIR,
        "tables",
        "gene_identifier_mapping.csv"
      )
    )


    # ========================================================
    # 3. GO BP AND KEGG OVER-REPRESENTATION ANALYSIS
    # ========================================================

    enrichment_classes <- de_significant |>
      dplyr::distinct(
        .data$cell_type,
        .data$direction,
        .data$class
      ) |>
      dplyr::arrange(
        match(
          .data$cell_type,
          de_cell_type_levels
        ),
        match(
          .data$direction,
          c(
            up_case_label,
            up_control_label
          )
        )
      )


    enrichment_audit_list <- list()
    go_result_list <- list()
    kegg_result_list <- list()


    # --------------------------------------------------------
    # Download the species-specific KEGG membership table ONCE
    # and reuse it for every DEG class. Calling enrichKEGG()
    # separately for each class can repeatedly contact the remote
    # service and may leave only one biological class missing after
    # a transient timeout. clusterProfiler::enricher() performs the
    # same over-representation test using the cached TERM2GENE and
    # TERM2NAME tables, while preserving each cell-type-specific
    # tested-gene universe.
    #
    # 物种特异KEGG成员表只下载一次，并供所有DEG类别共同使用。
    # 若对每个类别分别调用enrichKEGG()，会反复访问远程服务，短暂
    # 超时可能造成仅某一个生物学类别缺失。enricher()使用缓存的
    # TERM2GENE和TERM2NAME执行相同的过度富集检验，同时仍为每个
    # cell type保留其自身edgeR受检基因背景。
    # --------------------------------------------------------

    kegg_ora_cache_file <- file.path(
      MDIR,
      "tables",
      "KEGG_ORA_gene_sets.csv"
    )


    kegg_ora_gene_sets <- NULL
    kegg_annotation_error <- ""


    if (file.exists(kegg_ora_cache_file)) {
      cached_kegg_ora <- data.table::fread(
        kegg_ora_cache_file,
        data.table = FALSE
      )


      required_kegg_ora_columns <- c(
        "species",
        "kegg_code",
        "pathway_id",
        "pathway_name",
        "ENTREZID"
      )


      if (
        all(required_kegg_ora_columns %in% colnames(cached_kegg_ora)) &&
        any(cached_kegg_ora$species == species_key)
      ) {
        kegg_ora_gene_sets <- cached_kegg_ora |>
          dplyr::filter(
            .data$species == species_key
          )
      }
    }


    if (is.null(kegg_ora_gene_sets) || nrow(kegg_ora_gene_sets) == 0L) {

      old_timeout <- getOption("timeout")
      options(timeout = max(3600, old_timeout))


      kegg_download <- NULL


      for (kegg_attempt in seq_len(2L)) {
        kegg_download <- tryCatch(
          clusterProfiler::download_KEGG(
            species = species_cfg$kegg_code,
            keggType = "KEGG",
            keyType = "ncbi-geneid"
          ),
          error = function(e) {
            kegg_annotation_error <<- conditionMessage(e)
            NULL
          }
        )


        if (!is.null(kegg_download)) {
          kegg_annotation_error <- ""
          break
        }
      }


      options(timeout = old_timeout)


      if (!is.null(kegg_download)) {
        kegg_membership <- as.data.frame(
          kegg_download$KEGGPATHID2EXTID,
          stringsAsFactors = FALSE
        )[, seq_len(2L), drop = FALSE]


        colnames(kegg_membership) <- c(
          "pathway_id",
          "ENTREZID"
        )


        kegg_pathway_names <- as.data.frame(
          kegg_download$KEGGPATHID2NAME,
          stringsAsFactors = FALSE
        )[, seq_len(2L), drop = FALSE]


        colnames(kegg_pathway_names) <- c(
          "pathway_id",
          "pathway_name"
        )


        kegg_ora_gene_sets <- kegg_membership |>
          dplyr::mutate(
            pathway_id = sub(
              "^path:",
              "",
              as.character(.data$pathway_id)
            ),
            ENTREZID = sub(
              "^[^:]+:",
              "",
              as.character(.data$ENTREZID)
            )
          ) |>
          dplyr::inner_join(
            kegg_pathway_names |>
              dplyr::mutate(
                pathway_id = sub(
                  "^path:",
                  "",
                  as.character(.data$pathway_id)
                )
              ),
            by = "pathway_id"
          ) |>
          dplyr::mutate(
            species = species_key,
            kegg_code = species_cfg$kegg_code,
            .before = 1
          ) |>
          dplyr::distinct()


        data.table::fwrite(
          kegg_ora_gene_sets,
          kegg_ora_cache_file
        )
      } else {
        warning(
          paste0(
            "Species-specific KEGG annotation could not be downloaded after ",
            "two attempts: ",
            kegg_annotation_error,
            ". GO analysis will continue; KEGG status will be recorded."
          ),
          call. = FALSE
        )
      }
    }


    kegg_term2gene <- if (
      is.null(kegg_ora_gene_sets)
    ) {
      NULL
    } else {
      kegg_ora_gene_sets |>
        dplyr::select(
          dplyr::all_of(
            c(
              "pathway_id",
              "ENTREZID"
            )
          )
        ) |>
        dplyr::distinct()
    }


    kegg_term2name <- if (
      is.null(kegg_ora_gene_sets)
    ) {
      NULL
    } else {
      kegg_ora_gene_sets |>
        dplyr::select(
          dplyr::all_of(
            c(
              "pathway_id",
              "pathway_name"
            )
          )
        ) |>
        dplyr::distinct()
    }


    for (
      enrichment_index in
        seq_len(
          nrow(
            enrichment_classes
          )
        )
    ) {

      target_cell_type <- enrichment_classes$cell_type[[
        enrichment_index
      ]]


      target_direction <- enrichment_classes$direction[[
        enrichment_index
      ]]


      target_class <- enrichment_classes$class[[
        enrichment_index
      ]]


      target_symbols <- de_significant |>
        dplyr::filter(
          .data$cell_type == target_cell_type,
          .data$direction == target_direction
        ) |>
        dplyr::pull(
          .data$gene
        ) |>
        unique()


      universe_symbols <- de_all |>
        dplyr::filter(
          as.character(
            .data$cell_type
          ) == target_cell_type
        ) |>
        dplyr::pull(
          .data$gene
        ) |>
        unique()


      target_entrez <- identifier_mapping |>
        dplyr::filter(
          .data$SYMBOL %in%
            target_symbols
        ) |>
        dplyr::pull(
          .data$ENTREZID
        ) |>
        unique() |>
        as.character()


      universe_entrez <- identifier_mapping |>
        dplyr::filter(
          .data$SYMBOL %in%
            universe_symbols
        ) |>
        dplyr::pull(
          .data$ENTREZID
        ) |>
        unique() |>
        as.character()


      enrichment_status <- if (
        length(
          target_entrez
        ) <
          ENRICHMENT_MIN_MAPPED_DEGS
      ) {
        "SKIPPED"
      } else {
        "TESTED"
      }


      enrichment_reason <- if (
        enrichment_status == "SKIPPED"
      ) {
        paste0(
          "mapped significant genes ",
          length(
            target_entrez
          ),
          " < ",
          ENRICHMENT_MIN_MAPPED_DEGS
        )
      } else {
        ""
      }


      go_n_significant <- 0L
      kegg_n_significant <- 0L
      go_error_message <- ""
      kegg_error_message <- ""


      if (
        enrichment_status == "TESTED"
      ) {
        go_result <- tryCatch(
          clusterProfiler::enrichGO(
            gene = target_entrez,
            universe = universe_entrez,
            OrgDb = orgdb_object,
            keyType = "ENTREZID",
            ont = "BP",
            pAdjustMethod = "BH",
            pvalueCutoff = 1,
            qvalueCutoff = 1,
            readable = TRUE
          ),
          error = function(e) {
            go_error_message <<- conditionMessage(e)
            warning(
              paste0(
                "GO enrichment failed for ",
                target_class,
                ": ",
                conditionMessage(e)
              ),
              call. = FALSE
            )
            NULL
          }
        )


        if (
          !is.null(
            go_result
          )
        ) {
          go_table <- as.data.frame(
            go_result
          )


          if (
            nrow(
              go_table
            ) > 0L
          ) {
            go_table <- go_table |>
              dplyr::mutate(
                cell_type = target_cell_type,
                direction = target_direction,
                class = target_class,
                ontology = "GO Biological Process",
                .before = 1
              )

            go_result_list[[target_class]] <- go_table
            go_n_significant <- sum(
              go_table$p.adjust <=
                ENRICHMENT_FDR_THRESHOLD,
              na.rm = TRUE
            )
          }
        }


        kegg_result <- if (
          is.null(kegg_term2gene) ||
          is.null(kegg_term2name)
        ) {
          kegg_error_message <- paste0(
            "Species-specific KEGG annotation unavailable: ",
            kegg_annotation_error
          )
          NULL
        } else {
          kegg_result <- tryCatch(
            clusterProfiler::enricher(
              gene = target_entrez,
              universe = universe_entrez,
              pAdjustMethod = "BH",
              pvalueCutoff = 1,
              qvalueCutoff = 1,
              TERM2GENE = kegg_term2gene,
              TERM2NAME = kegg_term2name
            ),
            error = function(e) {
              kegg_error_message <<- conditionMessage(e)
              NULL
            }
          )
          kegg_result
        }


        if (
          is.null(
            kegg_result
          ) &&
          nzchar(
            kegg_error_message
          )
        ) {
          warning(
            paste0(
              "KEGG enrichment failed for ",
              target_class,
              ": ",
              kegg_error_message
            ),
            call. = FALSE
          )
        }


        if (
          !is.null(
            kegg_result
          )
        ) {
          kegg_table <- as.data.frame(
            kegg_result
          )


          if (
            nrow(
              kegg_table
            ) > 0L
          ) {
            kegg_table <- kegg_table |>
              dplyr::mutate(
                cell_type = target_cell_type,
                direction = target_direction,
                class = target_class,
                ontology = "KEGG Pathway",
                .before = 1
              )

            kegg_result_list[[target_class]] <- kegg_table
            kegg_n_significant <- sum(
              kegg_table$p.adjust <=
                ENRICHMENT_FDR_THRESHOLD,
              na.rm = TRUE
            )
          }
        }
      }


      enrichment_audit_list[[target_class]] <- data.frame(
        cell_type = target_cell_type,
        direction = target_direction,
        class = target_class,
        n_significant_symbols = length(
          target_symbols
        ),
        n_mapped_significant_entrez = length(
          target_entrez
        ),
        n_universe_symbols = length(
          universe_symbols
        ),
        n_mapped_universe_entrez = length(
          universe_entrez
        ),
        status = enrichment_status,
        reason = enrichment_reason,
        GO_BP_error = go_error_message,
        KEGG_error = kegg_error_message,
        n_significant_GO_BP_terms = go_n_significant,
        n_significant_KEGG_terms = kegg_n_significant,
        stringsAsFactors = FALSE
      )
    }


    enrichment_audit <- dplyr::bind_rows(
      enrichment_audit_list
    )


    go_all <- dplyr::bind_rows(
      go_result_list
    )


    kegg_all <- dplyr::bind_rows(
      kegg_result_list
    )


    data.table::fwrite(
      enrichment_audit,
      file.path(
        MDIR,
        "tables",
        "enrichment_class_audit.csv"
      )
    )


    data.table::fwrite(
      go_all,
      file.path(
        MDIR,
        "tables",
        "GO_BP_enrichment_all.csv"
      )
    )


    data.table::fwrite(
      kegg_all,
      file.path(
        MDIR,
        "tables",
        "KEGG_enrichment_all.csv"
      )
    )


    # ========================================================
    # 4. CELL-STYLE ENRICHMENT DOT PLOTS
    # ========================================================
    # The compact dot plot uses enrichment ratio on x, term on y, gene count as
    # size and -log10(FDR) as colour. Facets preserve the originating cell type
    # and direction instead of merging biologically distinct DEG classes.
    # 紧凑气泡图以富集比例为x轴、术语为y轴、基因数为气泡大小、
    # -log10(FDR)为颜色。分面保留原始cell type和变化方向，不合并生物学上
    # 不同的DEG类别。

    ratio_to_numeric <- function(x) {
      vapply(
        strsplit(
          as.character(x),
          "/",
          fixed = TRUE
        ),
        function(parts) {
          if (
            length(parts) != 2L
          ) {
            return(
              NA_real_
            )
          }
          as.numeric(parts[[1]]) /
            as.numeric(parts[[2]])
        },
        FUN.VALUE = numeric(1)
      )
    }


    save_enrichment_dotplot <- function(
        enrichment_table,
        ontology_name,
        file_name
    ) {
      if (
        nrow(
          enrichment_table
        ) == 0L
      ) {
        return(
          invisible(FALSE)
        )
      }


      plot_data <- enrichment_table |>
        dplyr::filter(
          is.finite(
            .data$p.adjust
          ),
          .data$p.adjust <=
            ENRICHMENT_FDR_THRESHOLD
        ) |>
        dplyr::group_by(
          .data$class
        ) |>
        dplyr::arrange(
          .data$p.adjust,
          dplyr::desc(
            .data$Count
          ),
          .by_group = TRUE
        ) |>
        dplyr::slice_head(
          n = as.integer(
            ENRICHMENT_TOP_TERMS_PER_CLASS
          )
        ) |>
        dplyr::ungroup()


      if (
        nrow(
          plot_data
        ) == 0L
      ) {
        return(
          invisible(FALSE)
        )
      }


      plot_data <- plot_data |>
        dplyr::mutate(
          enrichment_ratio = ratio_to_numeric(
            .data$GeneRatio
          ),
          negative_log10_FDR = -log10(
            pmax(
              .data$p.adjust,
              .Machine$double.xmin
            )
          ),
          display_term = paste0(
            .data$Description,
            "   "
          ),
          display_term = factor(
            .data$display_term,
            levels = rev(
              unique(
                .data$display_term
              )
            )
          ),
          class = factor(
            .data$class,
            levels = unique(
              enrichment_classes$class[
                enrichment_classes$class %in%
                  .data$class
              ]
            )
          )
        )


      enrichment_plot <- ggplot2::ggplot(
        plot_data,
        ggplot2::aes(
          x = .data$enrichment_ratio,
          y = .data$display_term,
          size = .data$Count,
          colour = .data$negative_log10_FDR
        )
      ) +
        ggplot2::geom_point(
          alpha = 0.90
        ) +
        ggplot2::scale_colour_gradient(
          low = "#4DBBD5",
          high = "#E64B35",
          name = expression(
            -log[10](FDR)
          )
        ) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(
            class
          ),
          scales = "free_y",
          space = "free_y",
          switch = "y"
        ) +
        ggplot2::labs(
          title = paste0(
            ontology_name,
            " enrichment of pseudobulk DEG classes"
          ),
          subtitle = paste0(
            "FDR <= ",
            ENRICHMENT_FDR_THRESHOLD,
            "; universe = genes tested by edgeR within each cell type"
          ),
          x = "Gene ratio",
          y = NULL,
          size = "DEG count"
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          strip.placement = "outside",
          strip.text.y.left = ggplot2::element_text(
            angle = 0,
            face = "bold"
          ),
          axis.text.y = ggplot2::element_text(
            size = 8
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          ),
          legend.position = "right"
        )


      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          file_name
        ),
        enrichment_plot,
        width = 13,
        height = max(
          7,
          0.30 * nrow(
            plot_data
          ) + 2.5
        ),
        dpi = 300
      )


      invisible(TRUE)
    }


    save_enrichment_dotplot(
      go_all,
      "GO Biological Process",
      "GO_BP_enrichment_dotplot.png"
    )


    save_enrichment_dotplot(
      kegg_all,
      "KEGG pathway",
      "KEGG_enrichment_dotplot.png"
    )


    # ========================================================
    # 5. PUBLICATION-STYLE ggkegg PATHWAY PANELS
    #    期刊级ggkegg通路映射组合图
    # ========================================================
    #
    # The enrichment test and pathway visualization answer different questions.
    # Enrichment uses the significant DEG set and an edgeR-tested universe. The
    # map below then highlights ONLY the significant DEGs contributing to that
    # enriched pathway. Non-significant tested genes are not coloured as if they
    # were discoveries. The original KEGG raster is overlaid so arrows, compound
    # labels, group boxes and graphical connections absent from the parsed edge
    # table remain visible.
    #
    # 富集检验与通路作图回答不同问题。富集使用显著DEG及当前cell type内
    # edgeR已检验基因背景；下方通路图只高亮实际贡献于该富集通路的
    # 显著DEG，不会将非显著基因着色得像新发现。同时叠加KEGG原始底图，
    # 保留只存在于原始图像或特殊KGML元素中的箭头、化合物、group及连线。

    ggkegg_audit_list <- list()
    ggkegg_mapped_deg_list <- list()


    if (
      RUN_GGKEGG &&
      nrow(
        kegg_all
      ) > 0L
    ) {
      ggkegg_directory <- file.path(
        MDIR,
        "figures",
        "ggkegg"
      )


      ggkegg_map_directory <- file.path(
        ggkegg_directory,
        "pathway_maps"
      )


      ggkegg_composite_directory <- file.path(
        ggkegg_directory,
        "composite_panels"
      )


      dir.create(
        ggkegg_map_directory,
        recursive = TRUE,
        showWarnings = FALSE
      )


      dir.create(
        ggkegg_composite_directory,
        recursive = TRUE,
        showWarnings = FALSE
      )


      # The selected pathways may change after data or thresholds change. Remove
      # only PNG files previously generated inside this dedicated module folder;
      # no table, checkpoint or output from another module is touched.
      # 数据或阈值改变后选中通路可能变化。重跑时仅清理该模块专用目录中
      # 旧PNG；不删除表格、checkpoint或其他模块的任何结果。
      previous_ggkegg_png <- list.files(
        ggkegg_directory,
        pattern = "\\.png$",
        recursive = TRUE,
        full.names = TRUE
      )


      if (
        length(
          previous_ggkegg_png
        ) > 0L
      ) {
        unlink(
          previous_ggkegg_png,
          force = TRUE
        )
      }


      # Select top pathways WITHIN each biological DEG class before applying the
      # global cap. This preserves representation across cell types and directions.
      # 先在每个生物学DEG类别内选取top pathway，再应用总数上限，从而保留
      # 不同cell type及变化方向在最终图片中的代表性。
      ggkegg_selection <- kegg_all |>
        dplyr::filter(
          is.finite(
            .data$p.adjust
          ),
          .data$p.adjust <=
            ENRICHMENT_FDR_THRESHOLD
        ) |>
        dplyr::group_by(
          .data$class
        ) |>
        dplyr::arrange(
          .data$p.adjust,
          dplyr::desc(
            .data$Count
          ),
          .by_group = TRUE
        ) |>
        dplyr::slice_head(
          n = as.integer(
            GGKEGG_PATHWAYS_PER_CLASS
          )
        ) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          .data$p.adjust,
          dplyr::desc(
            .data$Count
          )
        ) |>
        dplyr::slice_head(
          n = as.integer(
            GGKEGG_MAX_PATHWAYS
          )
        )


      data.table::fwrite(
        ggkegg_selection,
        file.path(
          MDIR,
          "tables",
          "ggkegg_pathway_selection.csv"
        )
      )


      ggkegg_logfc_limit <- max(
        1,
        min(
          4,
          stats::quantile(
            abs(
              de_significant$logFC
            ),
            probs = 0.98,
            na.rm = TRUE,
            names = FALSE
          )
        )
      )


      # Prepare one reusable UMAP table. Failure to find a UMAP does not remove
      # the pathway map; it only omits the optional location panel.
      # 一次性准备可重用UMAP数据。如果checkpoint中没有UMAP，通路图仍会
      # 正常生成，只会省略可选的cell type空间定位panel。
      ggkegg_umap_data <- tryCatch(
        {
          ggkegg_object <- load_checkpoint(
            "M11_annotated"
          )


          ggkegg_umap_name <- grep(
            "umap",
            names(
              ggkegg_object@reductions
            ),
            ignore.case = TRUE,
            value = TRUE
          )[[1]]


          ggkegg_embedding <- Seurat::Embeddings(
            ggkegg_object,
            reduction = ggkegg_umap_name
          )


          data.frame(
            UMAP_1 = ggkegg_embedding[
              ,
              1
            ],
            UMAP_2 = ggkegg_embedding[
              ,
              2
            ],
            cell_type = as.character(
              ggkegg_object$cell_type
            ),
            condition = as.character(
              ggkegg_object[[
                CONDITION_COL,
                drop = TRUE
              ]]
            ),
            stringsAsFactors = FALSE
          )
        },
        error = function(e) {
          NULL
        }
      )


      condition_colours <- c(
        stats::setNames(
          "#3C5488",
          CONTROL_GROUP
        ),
        stats::setNames(
          "#E64B35",
          CASE_GROUP
        )
      )


      for (
        ggkegg_index in
          seq_len(
            nrow(
              ggkegg_selection
            )
          )
      ) {
        pathway_row <- ggkegg_selection[
          ggkegg_index,
          ,
          drop = FALSE
        ]


        pathway_id <- as.character(
          pathway_row$ID[[1]]
        )


        pathway_cell_type <- as.character(
          pathway_row$cell_type[[1]]
        )


        pathway_direction <- as.character(
          pathway_row$direction[[1]]
        )


        pathway_description <- as.character(
          pathway_row$Description[[1]]
        )


        safe_pathway_name <- gsub(
          "[^A-Za-z0-9_-]+",
          "_",
          paste(
            pathway_cell_type,
            pathway_direction,
            pathway_id,
            sep = "_"
          )
        )


        pathway_status <- "OK"
        pathway_reason <- ""
        pathway_n_nodes <- NA_integer_
        pathway_n_edges <- NA_integer_
        pathway_n_isolated_nodes <- NA_integer_
        pathway_n_mapped_nodes <- NA_integer_
        pathway_n_mapped_degs <- 0L


        pathway_plot_result <- tryCatch(
          {
            pathway_entrez_ids <- unique(
              strsplit(
                as.character(
                  pathway_row$geneID[[1]]
                ),
                "/",
                fixed = TRUE
              )[[1]]
            )


            pathway_degs <- de_significant |>
              dplyr::filter(
                .data$cell_type == pathway_cell_type,
                .data$direction == pathway_direction
              ) |>
              dplyr::inner_join(
                identifier_mapping |>
                  dplyr::rename(
                    gene = "SYMBOL"
                  ),
                by = "gene"
              ) |>
              dplyr::filter(
                as.character(
                  .data$ENTREZID
                ) %in%
                  pathway_entrez_ids
              ) |>
              dplyr::group_by(
                .data$ENTREZID
              ) |>
              dplyr::arrange(
                .data$FDR,
                dplyr::desc(
                  abs(
                    .data$logFC
                  )
                ),
                .by_group = TRUE
              ) |>
              dplyr::slice_head(
                n = 1L
              ) |>
              dplyr::ungroup() |>
              dplyr::mutate(
                pathway_id = pathway_id,
                pathway_description = pathway_description,
                enrichment_FDR = as.numeric(
                  pathway_row$p.adjust[[1]]
                ),
                .before = 1
              )


            pathway_n_mapped_degs <- nrow(
              pathway_degs
            )


            if (
              pathway_n_mapped_degs == 0L
            ) {
              stop(
                paste0(
                  "No significant DEG could be mapped back to the enriched pathway / ",
                  "无显著DEG可映射回已富集通路。"
                )
              )
            }


            ggkegg_mapped_deg_list[[
              paste0(
                ggkegg_index,
                "_",
                pathway_id
              )
            ]] <- pathway_degs


            logfc_vector <- stats::setNames(
              as.numeric(
                pathway_degs$logFC
              ),
              paste0(
                species_cfg$kegg_code,
                ":",
                pathway_degs$ENTREZID
              )
            )


            pathway_graph <- ggkegg::pathway(
              pathway_id,
              use_cache = TRUE
            ) |>
              tidygraph::activate(
                "nodes"
              ) |>
              dplyr::mutate(
                DEG_logFC = ggkegg::node_numeric(
                  logfc_vector,
                  name = "name",
                  how = "any",
                  sep = " "
                ),
                is_mapped_DEG = is.finite(
                  .data$DEG_logFC
                )
              )


            pathway_node_table <- pathway_graph |>
              tidygraph::activate(
                "nodes"
              ) |>
              tibble::as_tibble()


            pathway_edge_table <- pathway_graph |>
              tidygraph::activate(
                "edges"
              ) |>
              tibble::as_tibble()


            pathway_degree <- igraph::degree(
              tidygraph::as.igraph(
                pathway_graph
              ),
              mode = "all"
            )


            pathway_n_nodes <- nrow(
              pathway_node_table
            )


            pathway_n_edges <- nrow(
              pathway_edge_table
            )


            pathway_n_isolated_nodes <- sum(
              pathway_degree == 0L
            )


            pathway_n_mapped_nodes <- sum(
              pathway_node_table$is_mapped_DEG,
              na.rm = TRUE
            )


            # The coloured rectangles are drawn first and the species-specific
            # KEGG image is overlaid afterwards. Transparent database background
            # colours reveal the DEG fill while preserving original text and lines.
            # 先绘制DEG着色矩形，再叠加物种特异KEGG原始图。将数据库底色透明化后，
            # 既能显示DEG颜色，又能保留原始基因标签、箭头和网络连线。
            pathway_plot <- ggraph::ggraph(
              pathway_graph,
              layout = "manual",
              x = x,
              y = y
            ) +
              ggkegg::geom_node_rect(
                ggplot2::aes(
                  fill = .data$DEG_logFC,
                  filter = .data$type == "gene" &
                    .data$is_mapped_DEG
                ),
                colour = "#1F1F1F",
                linewidth = 0.40
              ) +
              ggkegg::overlay_raw_map(
                pid = pathway_id,
                use_cache = TRUE,
                high_res = FALSE,
                transparent_colors = c(
                  "#FFFFFF",
                  "#BFBFFF",
                  "#BFFFBF"
                )
              ) +
              ggplot2::scale_fill_gradient2(
                low = "#3C5488",
                mid = "#F7F7F7",
                high = "#B2182B",
                midpoint = 0,
                limits = c(
                  -ggkegg_logfc_limit,
                  ggkegg_logfc_limit
                ),
                oob = scales::squish,
                name = paste0(
                  "DEG log2FC\n",
                  CASE_GROUP,
                  " vs ",
                  CONTROL_GROUP
                )
              ) +
              ggplot2::labs(
                title = pathway_description,
                subtitle = paste0(
                  pathway_cell_type,
                  " | ",
                  pathway_direction,
                  " | ",
                  pathway_id,
                  " | enrichment FDR = ",
                  formatC(
                    pathway_row$p.adjust[[1]],
                    format = "e",
                    digits = 2
                  )
                ),
                caption = paste0(
                  "Coloured nodes are significant pseudobulk DEG only; ",
                  "original KEGG labels and connections are retained."
                )
              ) +
              ggplot2::theme_void(
                base_size = 10
              ) +
              ggplot2::theme(
                plot.title = ggplot2::element_text(
                  face = "bold",
                  size = 12
                ),
                plot.subtitle = ggplot2::element_text(
                  size = 8.5,
                  colour = "#333333"
                ),
                plot.caption = ggplot2::element_text(
                  size = 7.5,
                  colour = "#4D4D4D",
                  hjust = 0
                ),
                legend.position = "right",
                plot.margin = ggplot2::margin(
                  10,
                  12,
                  8,
                  12
                )
              )


            save_plot_transparent(
              file.path(
                ggkegg_map_directory,
                paste0(
                  safe_pathway_name,
                  "_pathway.png"
                )
              ),
              pathway_plot,
              width = 13,
              height = 9,
              dpi = 300
            )


            pathway_degs_for_plot <- pathway_degs |>
              dplyr::arrange(
                .data$FDR,
                dplyr::desc(
                  abs(
                    .data$logFC
                  )
                )
              ) |>
              dplyr::slice_head(
                n = as.integer(
                  GGKEGG_MAX_GENES_IN_COMPANION
                )
              ) |>
              dplyr::mutate(
                gene = factor(
                  .data$gene,
                  levels = .data$gene[
                    order(
                      .data$logFC
                    )
                  ]
                ),
                negative_log10_FDR = -log10(
                  pmax(
                    .data$FDR,
                    .Machine$double.xmin
                  )
                )
              )


            gene_effect_plot <- ggplot2::ggplot(
              pathway_degs_for_plot,
              ggplot2::aes(
                x = .data$logFC,
                y = .data$gene
              )
            ) +
              ggplot2::geom_vline(
                xintercept = 0,
                colour = "#808080",
                linewidth = 0.35
              ) +
              ggplot2::geom_segment(
                ggplot2::aes(
                  x = 0,
                  xend = .data$logFC,
                  yend = .data$gene,
                  colour = .data$direction
                ),
                linewidth = 0.65,
                alpha = 0.75
              ) +
              ggplot2::geom_point(
                ggplot2::aes(
                  colour = .data$direction,
                  size = .data$negative_log10_FDR
                ),
                alpha = 0.95
              ) +
              ggplot2::scale_colour_manual(
                values = c(
                  stats::setNames(
                    "#B2182B",
                    up_case_label
                  ),
                  stats::setNames(
                    "#3C5488",
                    up_control_label
                  )
                ),
                drop = FALSE
              ) +
              ggplot2::scale_size_continuous(
                range = c(
                  2.0,
                  5.0
                ),
                name = expression(
                  -log[10](FDR)
                )
              ) +
              ggplot2::labs(
                title = "Mapped significant genes",
                x = paste0(
                  "log2FC (",
                  CASE_GROUP,
                  " vs ",
                  CONTROL_GROUP,
                  ")"
                ),
                y = NULL,
                colour = "Direction"
              ) +
              ggplot2::theme_classic(
                base_size = 9
              ) +
              ggplot2::theme(
                plot.title = ggplot2::element_text(
                  face = "bold",
                  size = 10
                ),
                axis.text.y = ggplot2::element_text(
                  face = "italic",
                  colour = "#222222"
                ),
                legend.position = "bottom",
                legend.box = "vertical"
              )


            pathway_umap_plot <- NULL


            if (
              !is.null(
                ggkegg_umap_data
              ) &&
              pathway_cell_type %in%
                ggkegg_umap_data$cell_type
            ) {
              target_umap_data <- ggkegg_umap_data |>
                dplyr::filter(
                  .data$cell_type == pathway_cell_type,
                  .data$condition %in%
                    c(
                      CONTROL_GROUP,
                      CASE_GROUP
                    )
                ) |>
                dplyr::mutate(
                  condition = factor(
                    .data$condition,
                    levels = c(
                      CONTROL_GROUP,
                      CASE_GROUP
                    )
                  )
                )


              pathway_umap_plot <- ggplot2::ggplot() +
                ggplot2::geom_point(
                  data = ggkegg_umap_data,
                  ggplot2::aes(
                    x = .data$UMAP_1,
                    y = .data$UMAP_2
                  ),
                  colour = "#D9D9D9",
                  size = 0.12,
                  alpha = 0.35
                ) +
                ggplot2::geom_point(
                  data = target_umap_data,
                  ggplot2::aes(
                    x = .data$UMAP_1,
                    y = .data$UMAP_2,
                    colour = .data$condition
                  ),
                  size = 0.32,
                  alpha = 0.85
                ) +
                ggplot2::scale_colour_manual(
                  values = condition_colours,
                  drop = FALSE
                ) +
                ggplot2::coord_equal() +
                ggplot2::labs(
                  title = pathway_cell_type,
                  subtitle = "Annotated cells in the global UMAP",
                  colour = "Condition"
                ) +
                ggplot2::theme_void(
                  base_size = 9
                ) +
                ggplot2::theme(
                  plot.title = ggplot2::element_text(
                    face = "bold",
                    size = 11
                  ),
                  plot.subtitle = ggplot2::element_text(
                    size = 8,
                    colour = "#4D4D4D"
                  ),
                  legend.position = "bottom"
                )
            }


            composite_plot <- if (
              is.null(
                pathway_umap_plot
              )
            ) {
              pathway_plot +
                gene_effect_plot +
                patchwork::plot_layout(
                  widths = c(
                    2.4,
                    1
                  ),
                  guides = "collect"
                )
            } else {
              pathway_umap_plot +
                pathway_plot +
                gene_effect_plot +
                patchwork::plot_layout(
                  widths = c(
                    0.85,
                    2.2,
                    0.95
                  ),
                  guides = "collect"
                )
            }


            composite_plot <- composite_plot &
              ggplot2::theme(
                legend.position = "bottom"
              )


            save_plot_transparent(
              file.path(
                ggkegg_composite_directory,
                paste0(
                  safe_pathway_name,
                  "_composite.png"
                )
              ),
              composite_plot,
              width = 18,
              height = 8.8,
              dpi = 300
            )


            TRUE
          },
          error = function(e) {
            pathway_status <<- "SKIPPED"
            pathway_reason <<- conditionMessage(e)
            FALSE
          }
        )


        ggkegg_audit_list[[
          paste0(
            ggkegg_index,
            "_",
            pathway_id
          )
        ]] <- data.frame(
          cell_type = pathway_cell_type,
          direction = pathway_direction,
          pathway_id = pathway_id,
          pathway_description = pathway_description,
          enrichment_FDR = pathway_row$p.adjust[[1]],
          n_KGML_nodes = pathway_n_nodes,
          n_KGML_edges = pathway_n_edges,
          n_isolated_KGML_nodes = pathway_n_isolated_nodes,
          n_significant_DEGs_mapped = pathway_n_mapped_degs,
          n_KEGG_nodes_with_mapped_DEG = pathway_n_mapped_nodes,
          original_KEGG_map_overlaid = pathway_status == "OK",
          status = pathway_status,
          reason = pathway_reason,
          stringsAsFactors = FALSE
        )
      }
    }


    ggkegg_audit <- if (
      length(
        ggkegg_audit_list
      ) == 0L
    ) {
      data.frame(
        cell_type = NA_character_,
        direction = NA_character_,
        pathway_id = NA_character_,
        pathway_description = NA_character_,
        enrichment_FDR = NA_real_,
        n_KGML_nodes = NA_integer_,
        n_KGML_edges = NA_integer_,
        n_isolated_KGML_nodes = NA_integer_,
        n_significant_DEGs_mapped = NA_integer_,
        n_KEGG_nodes_with_mapped_DEG = NA_integer_,
        original_KEGG_map_overlaid = FALSE,
        status = "SKIPPED",
        reason = if (
          !RUN_GGKEGG
        ) {
          "RUN_GGKEGG is FALSE"
        } else {
          "No FDR-significant KEGG pathway was available"
        },
        stringsAsFactors = FALSE
      )
    } else {
      dplyr::bind_rows(
        ggkegg_audit_list
      )
    }


    ggkegg_mapped_degs <- if (
      length(
        ggkegg_mapped_deg_list
      ) == 0L
    ) {
      data.frame(
        pathway_id = character(0),
        pathway_description = character(0),
        enrichment_FDR = numeric(0),
        cell_type = character(0),
        direction = character(0),
        gene = character(0),
        ENTREZID = character(0),
        logFC = numeric(0),
        FDR = numeric(0),
        stringsAsFactors = FALSE
      )
    } else {
      dplyr::bind_rows(
        ggkegg_mapped_deg_list
      )
    }


    data.table::fwrite(
      ggkegg_audit,
      file.path(
        MDIR,
        "tables",
        "ggkegg_pathway_audit.csv"
      )
    )


    data.table::fwrite(
      ggkegg_mapped_degs,
      file.path(
        MDIR,
        "tables",
        "ggkegg_mapped_significant_DEGs.csv"
      )
    )


    data.table::fwrite(
      data.frame(
        status = "OK",
        n_DEG_classes = nrow(
          enrichment_classes
        ),
        n_GO_BP_FDR_significant = if (
          nrow(go_all) == 0L
        ) 0L else sum(
          go_all$p.adjust <=
            ENRICHMENT_FDR_THRESHOLD,
          na.rm = TRUE
        ),
        n_KEGG_FDR_significant = if (
          nrow(kegg_all) == 0L
        ) 0L else sum(
          kegg_all$p.adjust <=
            ENRICHMENT_FDR_THRESHOLD,
          na.rm = TRUE
        ),
        interpretation = paste0(
          "Over-representation and pathway overlays are functional ",
          "interpretations of M13 pseudobulk DEG, not direct pathway-activity ",
          "or causal-signalling tests."
        ),
        stringsAsFactors = FALSE
      ),
      file.path(
        MDIR,
        "tables",
        "functional_interpretation_status.csv"
      )
    )


    write_module_status(
      MODULE,
      "OK"
    )
  }
}



# ============================================================
# M13C. REPLICATE-AWARE PSEUDOBULK GSVA
# M13C. 生物学重复感知的Pseudobulk GSVA
# ============================================================
#
# 【Module scope / 模块范围】
# This module estimates KEGG pathway-level expression scores separately for
# each annotated cell type. One biological replicate is one analysis column;
# individual cells are never treated as independent experimental replicates.
# 本模块在每个已注释细胞类型内分别估计KEGG通路层面的表达分数。每个生物学重复
# 对应一个分析列，绝不把单个细胞错误地当作彼此独立的实验重复。
#
# 【Input inherited from M13 / 继承自M13的输入】
# M13 saves the exact TMM-normalized logCPM matrix and aligned replicate
# metadata used after its cell-number, library-size and design checks. Reusing
# these objects keeps the DEG and GSVA analyses on the same eligible samples.
# M13保存通过细胞数、文库量和模型设计检查后的TMM标准化logCPM矩阵及其对齐的重复
# 信息。本模块直接复用这些对象，使差异基因分析与GSVA使用完全一致的合格样本。
#
# 【Statistical design / 统计设计】
# GSVA is applied to continuous logCPM values with a Gaussian kernel. Pathway
# scores are then tested by limma using the same condition contrast and optional
# covariates as M13. FDR correction is performed across all tested pathways
# within each cell type.
# GSVA对连续型logCPM值使用Gaussian核；随后使用limma以及与M13相同的条件比较和
# 可选协变量检验通路分数。每个细胞类型内部对全部受检通路进行FDR校正。
#
# 【Interpretation boundary / 解读边界】
# A GSVA score summarizes coordinated expression of genes in a pathway. It is
# not direct evidence of pathway flux, protein activity, metabolite abundance,
# receptor-ligand signalling or causality. Significant pathways should be read
# together with the mapped DEG, KEGG pathway maps and the underlying biology.
# GSVA分数反映通路基因的协调表达，不等同于真实代谢通量、蛋白活性、代谢物丰度、
# 受体-配体信号或因果关系。显著通路应结合映射到通路上的差异基因、KEGG通路图和
# 已知生物学共同解释。
#
# 【Principal output / 主要输出】
#
# tables/GSVA_all_cell_types_limma.csv
#   Complete pathway-level limma results / 完整通路层面limma结果
#
# tables/GSVA_significant_pathways.csv
#   FDR- and effect-size-filtered pathways / 经FDR及效应量筛选的通路
#
# tables/GSVA_scores_long.csv
#   Replicate-level pathway scores for reproducibility / 重复层面通路分数长表
#
# figures/GSVA_differential_pathways_dotplot.png
#   Cross-cell-type pathway effect overview / 跨细胞类型通路效应总览
#
# figures/GSVA_cross_cell_type_effect_heatmap.png
#   Heatmap of pathway score differences / 通路分数差异热图
#
# figures/sample_scores/ and figures/heatmaps/
#   Cell-type-specific replicate plots and heatmaps / 各细胞类型重复点图及热图
#

if (RUN_M13C_GSVA) {

  MODULE <- "M13C_pseudobulk_GSVA"
  MDIR <- module_dir(MODULE)


  required_packages <- c(
    "GSVA",
    "limma",
    "BiocParallel",
    "AnnotationDbi",
    "clusterProfiler",
    species_cfg$orgdb_package,
    "ComplexHeatmap",
    "circlize"
  )


  missing_packages <- required_packages[
    !vapply(
      required_packages,
      requireNamespace,
      logical(1),
      quietly = TRUE
    )
  ]


  if (length(missing_packages) > 0L) {
    safe_stop_module(
      MODULE,
      paste0(
        "Missing packages / 缺少软件包: ",
        paste(missing_packages, collapse = ", "),
        ". Re-run M00 with INSTALL_PACKAGES = TRUE. / ",
        "请设置INSTALL_PACKAGES = TRUE并重新运行M00。"
      )
    )
  }


  dir.create(
    file.path(MDIR, "figures", "sample_scores"),
    recursive = TRUE,
    showWarnings = FALSE
  )


  dir.create(
    file.path(MDIR, "figures", "heatmaps"),
    recursive = TRUE,
    showWarnings = FALSE
  )


  m13_table_dir <- file.path(
    DIR_RESULTS,
    "modules",
    "M13_pseudobulk_DE",
    "tables"
  )


  expression_file <- file.path(
    m13_table_dir,
    "pseudobulk_logCPM_by_cell_type.rds"
  )


  metadata_file <- file.path(
    m13_table_dir,
    "pseudobulk_sample_metadata_by_cell_type.rds"
  )


  if (!file.exists(expression_file) || !file.exists(metadata_file)) {
    safe_stop_module(
      MODULE,
      paste0(
        "M13 GSVA input is missing. Re-run M13 before M13C. / ",
        "缺少M13生成的GSVA输入，请先重新运行M13。\n",
        expression_file,
        "\n",
        metadata_file
      )
    )
  }


  pseudobulk_logcpm_list <- readRDS(expression_file)
  pseudobulk_metadata_list <- readRDS(metadata_file)


  common_cell_types <- intersect(
    names(pseudobulk_logcpm_list),
    names(pseudobulk_metadata_list)
  )


  if (length(common_cell_types) == 0L) {
    safe_stop_module(
      MODULE,
      paste0(
        "No eligible cell type was shared by the M13 expression and metadata ",
        "objects. / M13表达矩阵与样本信息中没有共同的合格细胞类型。"
      )
    )
  }


  # ----------------------------------------------------------
  # 1. BUILD OR REUSE THE SPECIES-SPECIFIC KEGG GENE SETS
  # 1. 构建或复用物种特异的KEGG基因集
  # ----------------------------------------------------------
  #
  # The complete pathway membership table is cached inside this module. This
  # makes later reruns reproducible and avoids unnecessary repeated KEGG calls.
  # 完整的通路成员表缓存在本模块内，使后续重跑可复现，并避免不必要地重复访问KEGG。

  kegg_gene_set_file <- file.path(
    MDIR,
    "tables",
    "KEGG_GSVA_gene_sets.csv"
  )


  kegg_gene_set_table <- NULL


  if (file.exists(kegg_gene_set_file)) {
    cached_kegg_gene_sets <- data.table::fread(
      kegg_gene_set_file,
      data.table = FALSE
    )


    required_gene_set_columns <- c(
      "species",
      "kegg_code",
      "pathway_id",
      "pathway_name",
      "ENTREZID",
      "SYMBOL"
    )


    if (
      all(required_gene_set_columns %in% colnames(cached_kegg_gene_sets)) &&
      any(cached_kegg_gene_sets$species == species_key)
    ) {
      kegg_gene_set_table <- cached_kegg_gene_sets |>
        dplyr::filter(.data$species == species_key)
    }
  }


  if (is.null(kegg_gene_set_table) || nrow(kegg_gene_set_table) == 0L) {

    old_timeout <- getOption("timeout")
    options(timeout = max(3600, old_timeout))


    kegg_download <- NULL
    kegg_download_error <- ""


    for (attempt in seq_len(2L)) {
      kegg_download <- tryCatch(
        clusterProfiler::download_KEGG(
          species = species_cfg$kegg_code,
          keggType = "KEGG",
          keyType = "ncbi-geneid"
        ),
        error = function(e) {
          kegg_download_error <<- conditionMessage(e)
          NULL
        }
      )


      if (!is.null(kegg_download)) {
        break
      }
    }


    options(timeout = old_timeout)


    if (is.null(kegg_download)) {
      safe_stop_module(
        MODULE,
        paste0(
          "KEGG gene-set download failed after two attempts: ",
          kegg_download_error,
          " / KEGG基因集连续两次下载失败。请检查网络后重新运行M13C；",
          "成功下载后文件会被缓存，无需每次联网。"
        )
      )
    }


    kegg_membership <- as.data.frame(
      kegg_download$KEGGPATHID2EXTID,
      stringsAsFactors = FALSE
    )[, seq_len(2L), drop = FALSE]


    colnames(kegg_membership) <- c(
      "pathway_id",
      "ENTREZID"
    )


    kegg_pathway_names <- as.data.frame(
      kegg_download$KEGGPATHID2NAME,
      stringsAsFactors = FALSE
    )[, seq_len(2L), drop = FALSE]


    colnames(kegg_pathway_names) <- c(
      "pathway_id",
      "pathway_name"
    )


    kegg_membership <- kegg_membership |>
      dplyr::mutate(
        pathway_id = sub("^path:", "", as.character(.data$pathway_id)),
        ENTREZID = sub("^[A-Za-z]+:", "", as.character(.data$ENTREZID))
      )


    kegg_pathway_names <- kegg_pathway_names |>
      dplyr::mutate(
        pathway_id = sub("^path:", "", as.character(.data$pathway_id)),
        pathway_name = as.character(.data$pathway_name)
      )


    orgdb <- getExportedValue(
      species_cfg$orgdb_package,
      species_cfg$orgdb_package
    )


    entrez_to_symbol <- suppressMessages(
      AnnotationDbi::select(
        orgdb,
        keys = unique(kegg_membership$ENTREZID),
        keytype = "ENTREZID",
        columns = "SYMBOL"
      )
    ) |>
      dplyr::filter(
        !is.na(.data$SYMBOL),
        nzchar(.data$SYMBOL)
      ) |>
      dplyr::distinct(
        .data$ENTREZID,
        .data$SYMBOL
      )


    kegg_gene_set_table <- kegg_membership |>
      dplyr::inner_join(
        kegg_pathway_names,
        by = "pathway_id"
      ) |>
      dplyr::inner_join(
        entrez_to_symbol,
        by = "ENTREZID"
      ) |>
      dplyr::mutate(
        species = species_key,
        kegg_code = species_cfg$kegg_code,
        .before = 1
      ) |>
      dplyr::distinct()


    data.table::fwrite(
      kegg_gene_set_table,
      kegg_gene_set_file
    )
  }


  kegg_pathway_metadata <- kegg_gene_set_table |>
    dplyr::distinct(
      .data$pathway_id,
      .data$pathway_name
    )


  kegg_gene_sets <- split(
    kegg_gene_set_table$SYMBOL,
    kegg_gene_set_table$pathway_id
  ) |>
    lapply(unique)


  # ----------------------------------------------------------
  # 2. CELL-TYPE-SPECIFIC GSVA AND LIMMA
  # 2. 各细胞类型内的GSVA与limma检验
  # ----------------------------------------------------------

  gsva_result_list <- list()
  gsva_score_list <- list()
  gsva_status_list <- list()


  condition_colours <- stats::setNames(
    c("#3C5488", "#E64B35"),
    c(CONTROL_GROUP, CASE_GROUP)
  )


  safe_gsva_celltype_names <- make.unique(
    gsub(
      "[^A-Za-z0-9]+",
      "_",
      common_cell_types
    ),
    sep = "_"
  )


  names(safe_gsva_celltype_names) <- common_cell_types


  for (target_cell_type in common_cell_types) {

    safe_cell_type <- safe_gsva_celltype_names[[target_cell_type]]
    expression_matrix <- as.matrix(
      pseudobulk_logcpm_list[[target_cell_type]]
    )


    replicate_metadata <- as.data.frame(
      pseudobulk_metadata_list[[target_cell_type]],
      stringsAsFactors = FALSE
    )


    if (
      is.null(rownames(expression_matrix)) ||
      is.null(colnames(expression_matrix)) ||
      nrow(expression_matrix) == 0L ||
      ncol(expression_matrix) == 0L
    ) {
      gsva_status_list[[target_cell_type]] <- data.frame(
        cell_type = target_cell_type,
        status = "SKIPPED",
        reason = "M13 logCPM matrix is empty or lacks dimnames",
        stringsAsFactors = FALSE
      )
      next
    }


    expression_matrix <- expression_matrix[
      !duplicated(rownames(expression_matrix)),
      ,
      drop = FALSE
    ]


    if (!"pseudobulk_id" %in% colnames(replicate_metadata)) {
      gsva_status_list[[target_cell_type]] <- data.frame(
        cell_type = target_cell_type,
        status = "SKIPPED",
        reason = "pseudobulk_id is absent from M13 metadata",
        stringsAsFactors = FALSE
      )
      next
    }


    replicate_metadata <- replicate_metadata[
      match(
        colnames(expression_matrix),
        replicate_metadata$pseudobulk_id
      ),
      ,
      drop = FALSE
    ]


    if (
      any(is.na(replicate_metadata$pseudobulk_id)) ||
      !identical(
        colnames(expression_matrix),
        as.character(replicate_metadata$pseudobulk_id)
      )
    ) {
      gsva_status_list[[target_cell_type]] <- data.frame(
        cell_type = target_cell_type,
        status = "SKIPPED",
        reason = "Expression columns and replicate metadata could not be aligned",
        stringsAsFactors = FALSE
      )
      next
    }


    replicate_metadata[[CONDITION_COL]] <- factor(
      as.character(replicate_metadata[[CONDITION_COL]]),
      levels = c(CONTROL_GROUP, CASE_GROUP)
    )


    model_variables <- unique(
      c(CONDITION_COL, DE_COVARIATES)
    )


    if (!all(model_variables %in% colnames(replicate_metadata))) {
      gsva_status_list[[target_cell_type]] <- data.frame(
        cell_type = target_cell_type,
        status = "SKIPPED",
        reason = paste0(
          "Missing model variables: ",
          paste(
            setdiff(model_variables, colnames(replicate_metadata)),
            collapse = ", "
          )
        ),
        stringsAsFactors = FALSE
      )
      next
    }


    design <- stats::model.matrix(
      stats::reformulate(
        model_variables,
        intercept = FALSE
      ),
      data = replicate_metadata
    )


    condition_control_column <- paste0(
      CONDITION_COL,
      CONTROL_GROUP
    )


    condition_case_column <- paste0(
      CONDITION_COL,
      CASE_GROUP
    )


    if (
      qr(design)$rank < ncol(design) ||
      nrow(design) <= ncol(design) ||
      !all(
        c(condition_control_column, condition_case_column) %in%
          colnames(design)
      )
    ) {
      gsva_status_list[[target_cell_type]] <- data.frame(
        cell_type = target_cell_type,
        status = "SKIPPED",
        reason = "The GSVA limma design is not estimable",
        stringsAsFactors = FALSE
      )
      next
    }


    pathway_gene_counts <- vapply(
      kegg_gene_sets,
      function(x) length(intersect(x, rownames(expression_matrix))),
      integer(1)
    )


    eligible_gene_sets <- kegg_gene_sets[
      pathway_gene_counts >= GSVA_MIN_GENE_SET_SIZE &
        pathway_gene_counts <= GSVA_MAX_GENE_SET_SIZE
    ]


    if (length(eligible_gene_sets) == 0L) {
      gsva_status_list[[target_cell_type]] <- data.frame(
        cell_type = target_cell_type,
        status = "SKIPPED",
        reason = "No KEGG gene set passed the configured size limits",
        stringsAsFactors = FALSE
      )
      next
    }


    gsva_parameter <- GSVA::gsvaParam(
      exprData = expression_matrix,
      geneSets = eligible_gene_sets,
      minSize = GSVA_MIN_GENE_SET_SIZE,
      maxSize = GSVA_MAX_GENE_SET_SIZE,
      kcdf = "Gaussian",
      maxDiff = TRUE,
      absRanking = FALSE,
      sparse = FALSE,
      checkNA = "yes"
    )


    pathway_score_matrix <- GSVA::gsva(
      gsva_parameter,
      verbose = FALSE,
      BPPARAM = BiocParallel::SerialParam(
        progressbar = FALSE
      )
    )


    if (inherits(pathway_score_matrix, "SummarizedExperiment")) {
      pathway_score_matrix <- SummarizedExperiment::assay(
        pathway_score_matrix
      )
    }


    pathway_score_matrix <- as.matrix(pathway_score_matrix)


    contrast_vector <- rep(
      0,
      ncol(design)
    )


    names(contrast_vector) <- colnames(design)
    contrast_vector[condition_case_column] <- 1
    contrast_vector[condition_control_column] <- -1


    contrast_name <- paste0(
      CASE_GROUP,
      "_vs_",
      CONTROL_GROUP
    )


    contrast_matrix <- matrix(
      contrast_vector,
      ncol = 1L,
      dimnames = list(
        names(contrast_vector),
        contrast_name
      )
    )


    gsva_fit <- limma::lmFit(
      pathway_score_matrix,
      design
    )


    gsva_fit <- limma::contrasts.fit(
      gsva_fit,
      contrasts = contrast_matrix
    )


    gsva_fit <- limma::eBayes(
      gsva_fit,
      trend = FALSE,
      robust = TRUE
    )


    gsva_result <- limma::topTable(
      gsva_fit,
      coef = 1L,
      number = Inf,
      sort.by = "P"
    ) |>
      tibble::rownames_to_column("pathway_id") |>
      dplyr::left_join(
        kegg_pathway_metadata,
        by = "pathway_id"
      ) |>
      dplyr::mutate(
        cell_type = target_cell_type,
        contrast = contrast_name,
        n_control = sum(
          replicate_metadata[[CONDITION_COL]] == CONTROL_GROUP,
          na.rm = TRUE
        ),
        n_case = sum(
          replicate_metadata[[CONDITION_COL]] == CASE_GROUP,
          na.rm = TRUE
        ),
        direction = dplyr::if_else(
          .data$logFC >= 0,
          paste0("Higher in ", CASE_GROUP),
          paste0("Higher in ", CONTROL_GROUP)
        ),
        significant = .data$adj.P.Val <= GSVA_FDR_THRESHOLD &
          abs(.data$logFC) >= GSVA_MIN_ABS_SCORE_DIFFERENCE,
        .before = 1
      )


    gsva_result_list[[target_cell_type]] <- gsva_result


    gsva_score_long <- as.data.frame(
      pathway_score_matrix,
      stringsAsFactors = FALSE
    ) |>
      tibble::rownames_to_column("pathway_id") |>
      tidyr::pivot_longer(
        cols = -"pathway_id",
        names_to = "pseudobulk_id",
        values_to = "GSVA_score"
      ) |>
      dplyr::left_join(
        replicate_metadata,
        by = "pseudobulk_id"
      ) |>
      dplyr::left_join(
        kegg_pathway_metadata,
        by = "pathway_id"
      ) |>
      dplyr::mutate(
        cell_type = target_cell_type,
        .before = 1
      )


    gsva_score_list[[target_cell_type]] <- gsva_score_long


    gsva_status_list[[target_cell_type]] <- data.frame(
      cell_type = target_cell_type,
      status = "OK",
      reason = "",
      n_genes = nrow(expression_matrix),
      n_replicates = ncol(expression_matrix),
      n_pathways_tested = nrow(gsva_result),
      n_significant_pathways = sum(gsva_result$significant, na.rm = TRUE),
      stringsAsFactors = FALSE
    )


    # Select pathways by FDR first and absolute effect second. This selection is
    # used only for visualization; the complete statistical table is retained.
    # 作图时优先按FDR、其次按绝对效应量筛选；完整统计结果仍全部保留在表格中。
    top_pathways <- gsva_result |>
      dplyr::arrange(
        .data$adj.P.Val,
        dplyr::desc(abs(.data$logFC))
      ) |>
      dplyr::slice_head(
        n = GSVA_TOP_PATHWAYS_PER_CELL_TYPE
      )


    top_score_data <- gsva_score_long |>
      dplyr::filter(
        .data$pathway_id %in% top_pathways$pathway_id
      ) |>
      dplyr::mutate(
        pathway_label = factor(
          .data$pathway_name,
          levels = rev(top_pathways$pathway_name)
        )
      )


    sample_score_plot <- ggplot2::ggplot(
      top_score_data,
      ggplot2::aes(
        x = .data[[CONDITION_COL]],
        y = .data$GSVA_score,
        colour = .data[[CONDITION_COL]]
      )
    ) +
      ggplot2::geom_boxplot(
        width = 0.56,
        outlier.shape = NA,
        linewidth = 0.45,
        alpha = 0.15
      ) +
      ggplot2::geom_point(
        position = ggplot2::position_jitter(
          width = 0.08,
          height = 0,
          seed = RANDOM_SEED
        ),
        size = 2.3,
        alpha = 0.9
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$pathway_label),
        scales = "free_y",
        ncol = 2
      ) +
      ggplot2::scale_colour_manual(
        values = condition_colours,
        drop = FALSE
      ) +
      ggplot2::labs(
        title = paste0(target_cell_type, ": replicate-level KEGG GSVA"),
        subtitle = paste0(
          "Each point is one biological replicate; pathways ranked by limma FDR"
        ),
        x = NULL,
        y = "GSVA score",
        colour = "Condition",
        caption = paste0(
          "Contrast: ",
          CASE_GROUP,
          " versus ",
          CONTROL_GROUP
        )
      ) +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(
        strip.text = ggplot2::element_text(
          face = "bold",
          size = 8.5
        ),
        legend.position = "top",
        axis.text.x = ggplot2::element_text(
          angle = 20,
          hjust = 1
        )
      )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "sample_scores",
        paste0(safe_cell_type, "_GSVA_sample_scores.png")
      ),
      sample_score_plot,
      width = 11,
      height = 3.1 + 2.2 * ceiling(nrow(top_pathways) / 2),
      units = "in",
      dpi = 300
    )


    top_score_matrix <- pathway_score_matrix[
      top_pathways$pathway_id,
      ,
      drop = FALSE
    ]


    display_names <- kegg_pathway_metadata$pathway_name[
      match(
        rownames(top_score_matrix),
        kegg_pathway_metadata$pathway_id
      )
    ]


    display_names[is.na(display_names)] <- rownames(top_score_matrix)[
      is.na(display_names)
    ]


    rownames(top_score_matrix) <- make.unique(display_names)


    top_score_z <- t(
      apply(
        top_score_matrix,
        1,
        function(x) {
          x_sd <- stats::sd(x, na.rm = TRUE)
          if (!is.finite(x_sd) || x_sd == 0) {
            return(rep(0, length(x)))
          }
          (x - mean(x, na.rm = TRUE)) / x_sd
        }
      )
    )


    top_score_z[!is.finite(top_score_z)] <- 0
    top_score_z <- pmax(pmin(top_score_z, 2), -2)


    column_condition <- as.character(
      replicate_metadata[[CONDITION_COL]]
    )


    gsva_column_annotation <- ComplexHeatmap::HeatmapAnnotation(
      Condition = column_condition,
      col = list(
        Condition = condition_colours
      ),
      annotation_name_gp = grid::gpar(
        fontface = "bold"
      )
    )


    gsva_heatmap <- ComplexHeatmap::Heatmap(
      top_score_z,
      name = "Row Z-score",
      col = circlize::colorRamp2(
        c(-2, 0, 2),
        c("#3C5488", "#F7F7F7", "#E64B35")
      ),
      cluster_rows = nrow(top_score_z) > 1L,
      cluster_columns = FALSE,
      show_row_names = TRUE,
      show_column_names = TRUE,
      column_labels = replicate_metadata$pseudobulk_id,
      column_split = factor(
        column_condition,
        levels = c(CONTROL_GROUP, CASE_GROUP)
      ),
      top_annotation = gsva_column_annotation,
      column_title = paste0(
        target_cell_type,
        ": replicate-level KEGG GSVA"
      ),
      row_names_gp = grid::gpar(fontsize = 8.5),
      column_names_gp = grid::gpar(fontsize = 8),
      border = TRUE,
      heatmap_legend_param = list(
        title = "Row Z-score"
      )
    )


    ragg::agg_png(
      filename = file.path(
        MDIR,
        "figures",
        "heatmaps",
        paste0(safe_cell_type, "_GSVA_sample_heatmap.png")
      ),
      width = 11,
      height = max(6, 0.48 * nrow(top_score_z) + 3.5),
      units = "in",
      res = 300,
      background = "transparent"
    )


    tryCatch(
      ComplexHeatmap::draw(
        gsva_heatmap,
        background = "transparent",
        heatmap_legend_side = "right",
        annotation_legend_side = "right"
      ),
      finally = {
        if (grDevices::dev.cur() > 1L) {
          grDevices::dev.off()
        }
      }
    )
  }


  gsva_status <- dplyr::bind_rows(gsva_status_list)


  data.table::fwrite(
    gsva_status,
    file.path(
      MDIR,
      "tables",
      "GSVA_status.csv"
    )
  )


  if (length(gsva_result_list) == 0L) {
    safe_stop_module(
      MODULE,
      paste0(
        "All cell types were skipped. Inspect GSVA_status.csv. / ",
        "所有细胞类型均被跳过，请查看GSVA_status.csv。"
      )
    )
  }


  gsva_all <- dplyr::bind_rows(gsva_result_list)
  gsva_significant <- gsva_all |>
    dplyr::filter(.data$significant)
  gsva_scores_all <- dplyr::bind_rows(gsva_score_list)


  data.table::fwrite(
    gsva_all,
    file.path(
      MDIR,
      "tables",
      "GSVA_all_cell_types_limma.csv"
    )
  )


  data.table::fwrite(
    gsva_significant,
    file.path(
      MDIR,
      "tables",
      "GSVA_significant_pathways.csv"
    )
  )


  data.table::fwrite(
    gsva_scores_all,
    file.path(
      MDIR,
      "tables",
      "GSVA_scores_long.csv"
    )
  )


  # ----------------------------------------------------------
  # 3. CROSS-CELL-TYPE PUBLICATION-STYLE SUMMARIES
  # 3. 跨细胞类型的出版级汇总图
  # ----------------------------------------------------------

  # Rank pathways globally by their best FDR across cell types, followed by the
  # largest absolute effect. The selected pathway set is then shown for every
  # cell type, producing a directly comparable matrix rather than disconnected
  # lists with different rows in each facet.
  # 首先按各通路在所有cell type中的最小FDR排序，再参考最大绝对效应量。
  # 随后在所有cell type中展示同一组通路，形成可直接横向比较的矩阵，而不是
  # 每个分面使用不同通路、难以对应的独立列表。
  gsva_overview_pathways <- gsva_all |>
    dplyr::group_by(
      .data$pathway_id,
      .data$pathway_name
    ) |>
    dplyr::summarise(
      minimum_FDR = min(.data$adj.P.Val, na.rm = TRUE),
      maximum_absolute_effect = max(abs(.data$logFC), na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      .data$minimum_FDR,
      dplyr::desc(.data$maximum_absolute_effect)
    ) |>
    dplyr::slice_head(
      n = GSVA_OVERVIEW_MAX_PATHWAYS
    ) |>
    dplyr::mutate(
      pathway_label = paste0(
        .data$pathway_name,
        " [",
        .data$pathway_id,
        "]"
      ),
      pathway_label = factor(
        .data$pathway_label,
        levels = rev(.data$pathway_label)
      )
    )


  selected_pathway_ids <- gsva_overview_pathways$pathway_id


  gsva_plot_selection <- gsva_all |>
    dplyr::filter(
      .data$pathway_id %in% selected_pathway_ids
    ) |>
    dplyr::mutate(
      minus_log10_FDR = -log10(
        pmax(.data$adj.P.Val, .Machine$double.xmin)
      ),
      pathway_label = factor(
        paste0(
          .data$pathway_name,
          " [",
          .data$pathway_id,
          "]"
        ),
        levels = levels(gsva_overview_pathways$pathway_label)
      ),
      cell_type = factor(
        .data$cell_type,
        levels = sort(unique(as.character(gsva_all$cell_type)))
      )
    )


  gsva_effect_limit <- max(
    abs(gsva_plot_selection$logFC),
    na.rm = TRUE
  )


  if (!is.finite(gsva_effect_limit) || gsva_effect_limit == 0) {
    gsva_effect_limit <- 1
  }


  gsva_dotplot <- ggplot2::ggplot(
    gsva_plot_selection,
    ggplot2::aes(
      x = .data$cell_type,
      y = .data$pathway_label,
      fill = .data$logFC,
      size = .data$minus_log10_FDR
    )
  ) +
    ggplot2::geom_point(
      shape = 21,
      colour = "#333333",
      stroke = 0.25,
      alpha = 0.95
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#3C5488",
      mid = "#F7F7F7",
      high = "#E64B35",
      midpoint = 0,
      limits = c(
        -gsva_effect_limit,
        gsva_effect_limit
      ),
      name = paste0(
        CASE_GROUP,
        " - ",
        CONTROL_GROUP
      )
    ) +
    ggplot2::scale_size_continuous(
      range = c(0.7, 6.5),
      name = expression(-log[10](FDR))
    ) +
    ggplot2::labs(
      title = "Differential KEGG pathway scores across cell types",
      subtitle = paste0(
        "Replicate-aware pseudobulk GSVA; colour shows effect and size shows FDR"
      ),
      x = NULL,
      y = NULL,
      caption = paste0(
        "Displayed: up to ",
        GSVA_OVERVIEW_MAX_PATHWAYS,
        " pathways ranked by minimum FDR and maximum absolute effect; ",
        "complete results are retained in the table"
      )
    ) +
    ggplot2::theme_minimal(base_size = 10.5) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_line(
        colour = "#E8E8E8",
        linewidth = 0.25
      ),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        face = "bold"
      ),
      legend.position = "right"
    )


  save_plot_transparent(
    file.path(
      MDIR,
      "figures",
      "GSVA_differential_pathways_dotplot.png"
    ),
    gsva_dotplot,
    width = max(
      13,
      0.62 * length(unique(gsva_plot_selection$cell_type)) + 7
    ),
    height = max(
      9,
      0.30 * length(selected_pathway_ids) + 3
    ),
    units = "in",
    dpi = 300,
    limitsize = FALSE
  )


  effect_heatmap_table <- gsva_all |>
    dplyr::filter(
      .data$pathway_id %in% selected_pathway_ids
    ) |>
    dplyr::select(
      dplyr::all_of(
        c(
          "pathway_id",
          "pathway_name",
          "cell_type",
          "logFC"
        )
      )
    ) |>
    tidyr::pivot_wider(
      names_from = "cell_type",
      values_from = "logFC"
    )


  effect_heatmap_matrix <- as.matrix(
    effect_heatmap_table[
      ,
      setdiff(
        colnames(effect_heatmap_table),
        c("pathway_id", "pathway_name")
      ),
      drop = FALSE
    ]
  )


  rownames(effect_heatmap_matrix) <- make.unique(
    paste0(
      effect_heatmap_table$pathway_name,
      " [",
      effect_heatmap_table$pathway_id,
      "]"
    )
  )


  maximum_effect <- max(
    abs(effect_heatmap_matrix),
    na.rm = TRUE
  )


  if (!is.finite(maximum_effect) || maximum_effect == 0) {
    maximum_effect <- 1
  }


  gsva_effect_heatmap <- ComplexHeatmap::Heatmap(
    effect_heatmap_matrix,
    name = paste0(CASE_GROUP, " - ", CONTROL_GROUP),
    col = circlize::colorRamp2(
      c(-maximum_effect, 0, maximum_effect),
      c("#3C5488", "#F7F7F7", "#E64B35")
    ),
    na_col = "#D9D9D9",
    cluster_rows = nrow(effect_heatmap_matrix) > 1L,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    show_column_names = TRUE,
    column_title = "Differential KEGG GSVA effect across cell types",
    row_names_gp = grid::gpar(fontsize = 7.5),
    column_names_gp = grid::gpar(
      fontsize = 9,
      fontface = "bold"
    ),
    border = TRUE,
    heatmap_legend_param = list(
      title = "GSVA effect"
    )
  )


  ragg::agg_png(
    filename = file.path(
      MDIR,
      "figures",
      "GSVA_cross_cell_type_effect_heatmap.png"
    ),
    width = max(10, 1.2 * ncol(effect_heatmap_matrix) + 7),
    height = max(9, 0.28 * nrow(effect_heatmap_matrix) + 4),
    units = "in",
    res = 300,
    background = "transparent"
  )


  tryCatch(
    ComplexHeatmap::draw(
      gsva_effect_heatmap,
      background = "transparent",
      heatmap_legend_side = "right"
    ),
    finally = {
      if (grDevices::dev.cur() > 1L) {
        grDevices::dev.off()
      }
    }
  )


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M14. OPTIONAL FOCUSED BIOLOGICAL ANALYSIS / 可选重点生物学分析
# ============================================================
#
# 【Module scope / 模块范围】
# This module prepares targeted visualizations or subsets only after global QC,
# clustering and annotation have been completed. It is intended for transparent
# follow-up of predefined or newly motivated biological questions.
# 本模块仅在全局QC、聚类和注释完成后生成重点可视化或子集，用于透明地跟进预先定义
# 或由全局结果提出的生物学问题。
#
# 【Interpretation boundary / 解读边界】
# Focused plots are descriptive unless linked to an appropriate replicate-level
# model. Selecting genes after viewing the same data is exploratory and should be
# labelled as such in reports; it does not constitute independent validation.
# 重点图形本身属于描述性结果，除非同时连接到合适的重复层面模型。若基因是在查看
# 同一数据后选择，则属于探索性分析，应在报告中明确说明，不能视为独立验证。
#
# 【Purpose / 目的】
#
# 完成无偏QC、聚类和annotation之后，
# 再针对具体研究假设查看指定genes或cell types。
#
#
# 【为什么放在最后？】
#
# 避免在预处理阶段因为预期结果而人为影响：
#
# QC
# clustering
# annotation
#
#
# 【Output / 输出】
#
# focus_gene_featureplots.png
# focus_cell_types_object.rds
#
# ============================================================


if (RUN_M14_FOCUSED_ANALYSIS) {

  MODULE <- "M14_focused_analysis"
  MDIR <- module_dir(MODULE)


  obj <- load_checkpoint(
    "M11_annotated"
  )


  ANALYSIS_ASSAY <- attr(
    obj,
    "ANALYSIS_ASSAY"
  )


  if (
    is.null(ANALYSIS_ASSAY) ||
    !ANALYSIS_ASSAY %in%
    SeuratObject::Assays(obj)
  ) {

    ANALYSIS_ASSAY <- DefaultAssay(
      obj
    )
  }


  DefaultAssay(
    obj
  ) <- ANALYSIS_ASSAY


  # ==========================================================
  # 1. FOCUS GENES
  # ==========================================================

  focus_genes_present <- intersect(
    FOCUS_GENES,
    rownames(obj)
  )


  focus_genes_missing <- setdiff(
    FOCUS_GENES,
    rownames(obj)
  )


  if (
    length(focus_genes_missing) > 0L
  ) {

    message(
      "Focus genes not found / 以下focus genes不存在:\n",
      paste(
        focus_genes_missing,
        collapse = ", "
      )
    )
  }


  if (
    length(focus_genes_present) > 0L
  ) {

    p_focus <- FeaturePlot(
      obj,
      features = focus_genes_present,
      reduction = "umap",
      ncol = 4,
      order = TRUE
    )


    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "focus_gene_featureplots.png"
      ),
      p_focus,
      width = 16,
      height = max(
        4,
        4 *
          ceiling(
            length(
              focus_genes_present
            ) / 4
          )
      ),
      dpi = 300
    )
  }


  # ==========================================================
  # 2. FOCUS CELL TYPES
  # ==========================================================

  if (
    length(
      FOCUS_CELL_TYPES
    ) > 0L
  ) {


    if (
      !"cell_type" %in%
      colnames(
        obj[[]]
      )
    ) {

      safe_stop_module(
        MODULE,
        "cell_type metadata is missing."
      )
    }


    focus_cells <- colnames(
      obj
    )[
      as.character(
        obj$cell_type
      ) %in%
        FOCUS_CELL_TYPES
    ]


    if (
      length(focus_cells) > 0L
    ) {

      focus_obj <- subset(
        obj,
        cells = focus_cells
      )


      saveRDS(
        focus_obj,
        file.path(
          MDIR,
          "focus_cell_types_object.rds"
        )
      )

    } else {

      message(
        "No requested focus cell types were found."
      )
    }
  }


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# M14B. OPTIONAL LINEAGE PSEUDOTIME / 可选谱系拟时分析
# ============================================================
#
# 【Module scope / 模块范围】
# This optional module reconstructs separate relative transcriptional trajectories
# for every biologically plausible profile listed in PSEUDOTIME_LINEAGE_PROFILES.
# It never pools all annotated cell types automatically. Each profile is subset,
# reprocessed, checked and reported independently so that unrelated lineages do
# not dominate within-lineage structure and one failed profile does not suppress
# valid output from the remaining profiles.
# 本可选模块会针对PSEUDOTIME_LINEAGE_PROFILES中每条具有合理生物学依据的谱系，
# 分别推断相对转录轨迹；绝不会自动合并全部已注释cell type。每条谱系均独立完成
# 子集提取、重新处理、质量判断和结果输出，既避免无关谱系掩盖谱系内部结构，也
# 避免某一条谱系失败时连带阻止其余有效谱系输出。
#
# 【Direction assumption / 方向假设】
# CONTROL_GROUP is used as the root-orientation anchor and CASE_GROUP as the
# terminal-orientation anchor. Root and terminal subclusters are selected from
# the expression-derived subset clusters, constrained by the root/endpoint cell
# identities and signatures configured for the current profile. Rank-based UCell
# signature evidence carries 50%, a compatible SingleR/ImmGen reference carries
# 25%, and condition enrichment carries 25%. When reference evidence is not
# applicable or unavailable, the remaining evidence is renormalized, retaining
# one-third condition influence. Condition labels orient the trajectory; they do
# not constitute observed time or independent evidence of disease progression.
# CONTROL_GROUP作为起点方向锚点，CASE_GROUP作为终点方向锚点。起点和终点从表达
# 数据形成的子集亚群中选择，并由当前profile配置的起点/终点cell type身份及标志
# 基因签名提供生物学约束。基于排名的UCell签名证据占50%，相容的SingleR/ImmGen
# 参考证据占25%，condition富集占25%。若参考证据不适用或不可用，其余证据会重新
# 归一化，condition仍保留三分之一影响。condition标签只用于确定推断方向，不是
# 真实采样时间，也不能作为疾病进展的独立证据。
#
# 【Method / 方法】
# 1. Subset the requested lineage and retain only CONTROL_GROUP/CASE_GROUP cells.
# 2. Recompute HVGs, PCA, neighbors, subclusters and UMAP inside the lineage.
# 3. Select sufficiently large root/terminal subclusters using biological identity
#    constraints plus weighted UCell, compatible reference-label and condition-
#    enrichment evidence.
# 4. Fit a Slingshot principal-curve trajectory in PCA space.
# 5. Report cell-level pseudotime, replicate summaries, marker trends and an
#    exploratory pseudotime-associated heatmap.
#
# 1. 提取指定谱系，仅保留CONTROL_GROUP与CASE_GROUP细胞。
# 2. 在谱系内部重新计算HVG、PCA、邻居图、亚群和UMAP。
# 3. 结合生物学身份约束，按设定权重综合UCell、相容参考标签及control/case富集
#    证据，选择足够大的起点/终点亚群。
# 4. 在PCA空间使用Slingshot主曲线拟合轨迹。
# 5. 输出cell层面拟时、重复层面汇总、marker趋势和探索性拟时相关热图。
#
# 【Interpretation boundary / 解读边界】
# Pseudotime is a relative ordering of transcriptional states. It cannot alone
# prove that one mature cell type differentiates into another, that control cells
# are chronologically earlier, or that treatment caused the inferred transition.
# Because condition contributes to root/terminal selection, downstream differences
# between conditions are descriptive and must not be presented as an independent
# hypothesis test.
# 拟时只是转录状态的相对顺序，不能单独证明一个成熟cell type分化为另一个cell
# type、control在真实时间上更早，或处理导致了推断转换。由于condition参与起点/
# 终点选择，后续组间差异属于描述性结果，不能作为独立假设检验汇报。
#
# 【Base summary output / 总汇总输出】
#   M14B_pseudotime/tables/pseudotime_profile_status.csv
#   M14B_pseudotime/tables/combined_pseudotime_direction_summary.csv
#
# 【Per-profile output / 每条谱系输出】
# Each profile uses the same P##_profile_id folder name under tables, figures and
# objects. This keeps the M14B module compact while preserving complete separation
# among the five trajectories.
# 每条谱系在tables、figures和objects下使用相同的P##_profile_id子目录名，在保持
# M14B结构紧凑的同时，确保五条轨迹的结果彼此完全分离。
#
#   tables/P##_profile_id/pseudotime_cell_metadata.csv
#   tables/P##_profile_id/lineage_subcluster_markers_all.csv
#   tables/P##_profile_id/lineage_subcluster_top_markers.csv
#   tables/P##_profile_id/lineage_subcluster_annotation_template.csv
#   tables/P##_profile_id/pseudotime_by_biological_replicate.csv
#   tables/P##_profile_id/pseudotime_cluster_selection_audit.csv
#   tables/P##_profile_id/pseudotime_endpoint_marker_audit.csv
#   tables/P##_profile_id/pseudotime_reference_label_audit.csv
#   tables/P##_profile_id/pseudotime_selection_evidence_audit.csv
#   tables/P##_profile_id/pseudotime_direction_summary.csv
#   tables/P##_profile_id/pseudotime_dynamic_gene_correlations.csv
#   tables/P##_profile_id/pseudotime_heatmap_matrix.csv
#   figures/P##_profile_id/pseudotime_trajectory_overview.png
#   figures/P##_profile_id/lineage_subcluster_marker_dotplot.png
#   figures/P##_profile_id/lineage_subcluster_marker_heatmap.png
#   figures/P##_profile_id/pseudotime_by_condition_and_replicate.png
#   figures/P##_profile_id/pseudotime_feature_gene_trends.png
#   figures/P##_profile_id/pseudotime_dynamic_gene_heatmap.png
#   objects/P##_profile_id/lineage_pseudotime_object.rds
# ============================================================


if (RUN_M14B_PSEUDOTIME) {

  BASE_PSEUDOTIME_MODULE <- "M14B_pseudotime"
  BASE_PSEUDOTIME_MDIR <- module_dir(
    BASE_PSEUDOTIME_MODULE
  )


  if (
    !is.list(
      PSEUDOTIME_LINEAGE_PROFILES
    ) ||
    length(
      PSEUDOTIME_LINEAGE_PROFILES
    ) == 0L
  ) {
    safe_stop_module(
      BASE_PSEUDOTIME_MODULE,
      paste0(
        "PSEUDOTIME_LINEAGE_PROFILES is empty / ",
        "PSEUDOTIME_LINEAGE_PROFILES为空。"
      )
    )
  }


  pseudotime_profile_status_list <- list()


  PSEUDOTIME_DEFAULT_N_HVG <- PSEUDOTIME_N_HVG
  PSEUDOTIME_DEFAULT_N_PCS <- PSEUDOTIME_N_PCS
  PSEUDOTIME_DEFAULT_CLUSTER_RESOLUTION <-
    PSEUDOTIME_CLUSTER_RESOLUTION
  PSEUDOTIME_DEFAULT_MIN_TOTAL_CELLS <-
    PSEUDOTIME_MIN_TOTAL_CELLS
  PSEUDOTIME_DEFAULT_MIN_CELLS_PER_CONDITION <-
    PSEUDOTIME_MIN_CELLS_PER_CONDITION
  PSEUDOTIME_DEFAULT_MIN_CLUSTER_CELLS <-
    PSEUDOTIME_MIN_CLUSTER_CELLS

  # The ImmGen reference is loaded at most once and reused by every compatible
  # profile. A persistent error message prevents repeated download attempts and
  # is written into each affected profile's audit table.
  # ImmGen参考在一次M14B运行中最多加载一次，并由所有相容profile共用。若加载
  # 失败，会保留错误信息，避免重复下载，并写入受影响profile的审计表。
  PSEUDOTIME_IMMGEN_REFERENCE <- NULL
  PSEUDOTIME_IMMGEN_REFERENCE_ERROR <- NA_character_


  for (
    pseudotime_profile_index in
      seq_along(
        PSEUDOTIME_LINEAGE_PROFILES
      )
  ) {

    pseudotime_profile <-
      PSEUDOTIME_LINEAGE_PROFILES[[
        pseudotime_profile_index
      ]]


    profile_value <- function(
        field,
        default
    ) {
      value <- pseudotime_profile[[field]]

      if (
        is.null(value)
      ) {
        default
      } else {
        value
      }
    }


    PSEUDOTIME_PROFILE_PRIORITY <- as.integer(
      profile_value(
        "priority",
        pseudotime_profile_index
      )
    )


    PSEUDOTIME_PROFILE_ID <- gsub(
      "[^A-Za-z0-9_]+",
      "_",
      as.character(
        profile_value(
          "profile_id",
          paste0(
            "lineage_",
            pseudotime_profile_index
          )
        )
      )
    )


    PSEUDOTIME_PROFILE_DISPLAY_NAME <- as.character(
      profile_value(
        "display_name",
        PSEUDOTIME_PROFILE_ID
      )
    )


    PSEUDOTIME_TARGET_CELL_TYPES <- unique(
      as.character(
        profile_value(
          "target_cell_types",
          character(0)
        )
      )
    )


    PSEUDOTIME_ROOT_CELL_TYPES <- unique(
      as.character(
        profile_value(
          "root_cell_types",
          character(0)
        )
      )
    )


    PSEUDOTIME_ENDPOINT_CELL_TYPES <- unique(
      as.character(
        profile_value(
          "endpoint_cell_types",
          character(0)
        )
      )
    )


    PSEUDOTIME_FEATURE_GENES <- unique(
      as.character(
        profile_value(
          "feature_genes",
          character(0)
        )
      )
    )


    PSEUDOTIME_ROOT_MARKER_GENES <- unique(
      as.character(
        profile_value(
          "root_marker_genes",
          character(0)
        )
      )
    )


    PSEUDOTIME_ENDPOINT_MARKER_GENES <- unique(
      as.character(
        profile_value(
          "endpoint_marker_genes",
          character(0)
        )
      )
    )


    PSEUDOTIME_REFERENCE_ROOT_LABEL_PATTERNS <- unique(
      as.character(
        profile_value(
          "reference_root_label_patterns",
          character(0)
        )
      )
    )


    PSEUDOTIME_REFERENCE_ENDPOINT_LABEL_PATTERNS <- unique(
      as.character(
        profile_value(
          "reference_endpoint_label_patterns",
          character(0)
        )
      )
    )


    PSEUDOTIME_ALLOW_EXPLORATORY_LOW_SUPPORT <- isTRUE(
      profile_value(
        "allow_exploratory_low_support",
        FALSE
      )
    )


    PSEUDOTIME_FORCE_EXPLORATORY <- isTRUE(
      profile_value(
        "force_exploratory",
        FALSE
      )
    )


    PSEUDOTIME_N_HVG_PROFILE <- as.integer(
      profile_value(
        "n_hvg",
        PSEUDOTIME_DEFAULT_N_HVG
      )
    )


    PSEUDOTIME_N_PCS_PROFILE <- as.integer(
      profile_value(
        "n_pcs",
        PSEUDOTIME_DEFAULT_N_PCS
      )
    )


    PSEUDOTIME_CLUSTER_RESOLUTION_PROFILE <- as.numeric(
      profile_value(
        "cluster_resolution",
        PSEUDOTIME_DEFAULT_CLUSTER_RESOLUTION
      )
    )


    PSEUDOTIME_MIN_TOTAL_CELLS_PROFILE <- as.integer(
      profile_value(
        "min_total_cells",
        PSEUDOTIME_DEFAULT_MIN_TOTAL_CELLS
      )
    )


    PSEUDOTIME_MIN_CELLS_PER_CONDITION_PROFILE <- as.integer(
      profile_value(
        "min_cells_per_condition",
        PSEUDOTIME_DEFAULT_MIN_CELLS_PER_CONDITION
      )
    )


    PSEUDOTIME_MIN_CLUSTER_CELLS_PROFILE <- as.integer(
      profile_value(
        "min_cluster_cells",
        PSEUDOTIME_DEFAULT_MIN_CLUSTER_CELLS
      )
    )


    PSEUDOTIME_N_HVG <- PSEUDOTIME_N_HVG_PROFILE
    PSEUDOTIME_N_PCS <- PSEUDOTIME_N_PCS_PROFILE
    PSEUDOTIME_CLUSTER_RESOLUTION <-
      PSEUDOTIME_CLUSTER_RESOLUTION_PROFILE
    PSEUDOTIME_MIN_TOTAL_CELLS <-
      PSEUDOTIME_MIN_TOTAL_CELLS_PROFILE
    PSEUDOTIME_MIN_CELLS_PER_CONDITION <-
      PSEUDOTIME_MIN_CELLS_PER_CONDITION_PROFILE
    PSEUDOTIME_MIN_CLUSTER_CELLS <-
      PSEUDOTIME_MIN_CLUSTER_CELLS_PROFILE


    PSEUDOTIME_PROFILE_DIRECTORY <- paste0(
      "P",
      sprintf(
        "%02d",
        PSEUDOTIME_PROFILE_PRIORITY
      ),
      "_",
      PSEUDOTIME_PROFILE_ID
    )


    # All lineage-specific files remain inside the single M14B module. Tables,
    # figures and serialized objects use matching profile subdirectories so that
    # outputs from five independent trajectories cannot overwrite one another.
    # 所有谱系文件均保存在同一个M14B模块内。表格、图片和序列化对象使用名称一致
    # 的profile子目录，避免五条独立轨迹的同名文件相互覆盖。
    MODULE <- BASE_PSEUDOTIME_MODULE
    MDIR <- BASE_PSEUDOTIME_MDIR
    PSEUDOTIME_TABLE_DIR <- file.path(
      MDIR,
      "tables",
      PSEUDOTIME_PROFILE_DIRECTORY
    )
    PSEUDOTIME_FIGURE_DIR <- file.path(
      MDIR,
      "figures",
      PSEUDOTIME_PROFILE_DIRECTORY
    )
    PSEUDOTIME_OBJECT_DIR <- file.path(
      MDIR,
      "objects",
      PSEUDOTIME_PROFILE_DIRECTORY
    )

    invisible(lapply(
      c(
        PSEUDOTIME_TABLE_DIR,
        PSEUDOTIME_FIGURE_DIR,
        PSEUDOTIME_OBJECT_DIR
      ),
      dir.create,
      recursive = TRUE,
      showWarnings = FALSE
    ))


    message(
      "Running pseudotime profile / 正在运行拟时谱系: ",
      PSEUDOTIME_PROFILE_DISPLAY_NAME
    )


    PSEUDOTIME_ANALYSIS_GRADE <- NA_character_
    lineage_cells <- character(0)


    pseudotime_profile_outcome <- tryCatch(
      {


  # ==========================================================
  # 1. PACKAGE AND CONFIGURATION VALIDATION
  # ==========================================================

  if (
    !requireNamespace(
      "slingshot",
      quietly = TRUE
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "The optional pseudotime module requires the Bioconductor package ",
        "slingshot / 可选拟时模块需要Bioconductor包slingshot。\n",
        "Run M00 with INSTALL_PACKAGES = TRUE, then rerun M14."
      )
    )
  }


  if (
    !requireNamespace(
      "UCell",
      quietly = TRUE
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "The pseudotime endpoint audit requires the Bioconductor package ",
        "UCell / 拟时起终点审计需要Bioconductor包UCell。\n",
        "Run M00 with INSTALL_PACKAGES = TRUE, then rerun M14."
      )
    )
  }


  if (
    length(
      PSEUDOTIME_TARGET_CELL_TYPES
    ) == 0L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "PSEUDOTIME_TARGET_CELL_TYPES is empty / ",
        "PSEUDOTIME_TARGET_CELL_TYPES为空。\n",
        "Specify one biologically plausible lineage; the module will not use ",
        "all cell types automatically / 请明确指定一个合理谱系，本模块不会自动使用",
        "全部cell type。"
      )
    )
  }


  if (
    !is.numeric(
      PSEUDOTIME_CLUSTER_RESOLUTION
    ) ||
    length(
      PSEUDOTIME_CLUSTER_RESOLUTION
    ) != 1L ||
    !is.finite(
      PSEUDOTIME_CLUSTER_RESOLUTION
    ) ||
    PSEUDOTIME_CLUSTER_RESOLUTION <= 0
  ) {
    safe_stop_module(
      MODULE,
      "PSEUDOTIME_CLUSTER_RESOLUTION must be one positive number."
    )
  }


  integer_settings <- c(
    PSEUDOTIME_N_HVG,
    PSEUDOTIME_N_PCS,
    PSEUDOTIME_MIN_TOTAL_CELLS,
    PSEUDOTIME_MIN_CELLS_PER_CONDITION,
    PSEUDOTIME_MIN_CLUSTER_CELLS,
    PSEUDOTIME_HEATMAP_TOP_GENES,
    PSEUDOTIME_HEATMAP_BINS
  )


  if (
    any(
      !is.finite(
        integer_settings
      )
    ) ||
    any(
      integer_settings < 2L
    )
  ) {
    safe_stop_module(
      MODULE,
      "All pseudotime integer settings must be finite and at least 2."
    )
  }


  if (
    !is.numeric(
      PSEUDOTIME_MIN_REPLICATE_COVERAGE
    ) ||
    length(
      PSEUDOTIME_MIN_REPLICATE_COVERAGE
    ) != 1L ||
    !is.finite(
      PSEUDOTIME_MIN_REPLICATE_COVERAGE
    ) ||
    PSEUDOTIME_MIN_REPLICATE_COVERAGE <= 0 ||
    PSEUDOTIME_MIN_REPLICATE_COVERAGE > 1
  ) {
    safe_stop_module(
      MODULE,
      "PSEUDOTIME_MIN_REPLICATE_COVERAGE must be in (0, 1]."
    )
  }


  if (
    !is.numeric(
      PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE
    ) ||
    length(
      PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE
    ) != 1L ||
    !is.finite(
      PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE
    ) ||
    PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE < 0 ||
    PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE > 1
  ) {
    safe_stop_module(
      MODULE,
      "PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE must be in [0, 1]."
    )
  }


  selection_weights <- c(
    signature = PSEUDOTIME_SELECTION_WEIGHT_SIGNATURE,
    reference = PSEUDOTIME_SELECTION_WEIGHT_REFERENCE,
    condition = PSEUDOTIME_SELECTION_WEIGHT_CONDITION
  )


  if (
    any(
      !is.finite(
        selection_weights
      )
    ) ||
    any(
      selection_weights < 0
    ) ||
    abs(
      sum(
        selection_weights
      ) - 1
    ) > 1e-8 ||
    selection_weights[["signature"]] <= 0 ||
    selection_weights[["condition"]] <= 0
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Pseudotime selection weights must be non-negative, sum to 1, and ",
        "retain positive signature and condition weights / 拟时选择权重必须为",
        "非负数、总和为1，且signature与condition权重必须大于0。"
      )
    )
  }


  if (
    !is.numeric(
      PSEUDOTIME_UCELL_MAX_RANK
    ) ||
    length(
      PSEUDOTIME_UCELL_MAX_RANK
    ) != 1L ||
    !is.finite(
      PSEUDOTIME_UCELL_MAX_RANK
    ) ||
    PSEUDOTIME_UCELL_MAX_RANK < 100L
  ) {
    safe_stop_module(
      MODULE,
      "PSEUDOTIME_UCELL_MAX_RANK must be one finite integer >= 100."
    )
  }


  # ==========================================================
  # 2. LOAD AND VALIDATE ANNOTATED OBJECT
  # ==========================================================

  obj <- load_checkpoint(
    "M11_annotated"
  )


  validate_complete_annotation(
    obj,
    MODULE
  )


  required_trajectory_meta <- c(
    BIOLOGICAL_REPLICATE_COL,
    CONDITION_COL,
    "cell_type"
  )


  missing_trajectory_meta <- setdiff(
    required_trajectory_meta,
    colnames(
      obj[[]]
    )
  )


  if (
    length(
      missing_trajectory_meta
    ) > 0L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Trajectory metadata is missing / 拟时所需metadata缺失: ",
        paste(
          missing_trajectory_meta,
          collapse = ", "
        )
      )
    )
  }


  observed_cell_types <- sort(
    unique(
      as.character(
        obj$cell_type
      )
    )
  )


  requested_target_types <-
    PSEUDOTIME_TARGET_CELL_TYPES


  missing_target_types <- setdiff(
    requested_target_types,
    observed_cell_types
  )


  PSEUDOTIME_TARGET_CELL_TYPES <- intersect(
    requested_target_types,
    observed_cell_types
  )


  if (
    length(
      missing_target_types
    ) > 0L
  ) {
    message(
      "Requested profile cell types not found and omitted / ",
      "以下预设cell type不存在，已略过: ",
      paste(
        missing_target_types,
        collapse = ", "
      )
    )
  }


  if (
    length(
      PSEUDOTIME_TARGET_CELL_TYPES
    ) == 0L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "None of the requested profile cell types were found / ",
        "该预设谱系的目标cell type均不存在。 Requested: ",
        paste(
          requested_target_types,
          collapse = ", "
        )
      )
    )
  }


  PSEUDOTIME_ROOT_CELL_TYPES <- intersect(
    PSEUDOTIME_ROOT_CELL_TYPES,
    PSEUDOTIME_TARGET_CELL_TYPES
  )


  PSEUDOTIME_ENDPOINT_CELL_TYPES <- intersect(
    PSEUDOTIME_ENDPOINT_CELL_TYPES,
    PSEUDOTIME_TARGET_CELL_TYPES
  )


  obj <- join_assay_layers_if_needed(
    obj,
    assay = "RNA"
  )


  if (
    !"RNA" %in%
    SeuratObject::Assays(obj)
  ) {
    safe_stop_module(
      MODULE,
      "RNA assay is required for pseudotime / 拟时分析需要RNA assay。"
    )
  }


  rna_layers <- SeuratObject::Layers(
    obj[["RNA"]]
  )


  if (
    !all(
      c(
        "counts",
        "data"
      ) %in%
        rna_layers
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "RNA counts and normalized data layers are both required / ",
        "拟时分析同时需要RNA counts与标准化data layer。"
      )
    )
  }


  # ==========================================================
  # 3. SELECT THE REQUESTED LINEAGE
  # ==========================================================

  obj_meta <- obj[[]]


  lineage_keep <-
    as.character(
      obj_meta$cell_type
    ) %in%
      PSEUDOTIME_TARGET_CELL_TYPES &
    as.character(
      obj_meta[[CONDITION_COL]]
    ) %in%
      c(
        CONTROL_GROUP,
        CASE_GROUP
      ) &
    !is.na(
      obj_meta[[BIOLOGICAL_REPLICATE_COL]]
    ) &
    trimws(
      as.character(
        obj_meta[[BIOLOGICAL_REPLICATE_COL]]
      )
    ) != ""


  lineage_cells <- rownames(
    obj_meta
  )[
    lineage_keep
  ]


  if (
    length(
      lineage_cells
    ) <
    PSEUDOTIME_MIN_TOTAL_CELLS
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Too few target-lineage cells / 目标谱系细胞数不足: ",
        length(lineage_cells),
        " < ",
        PSEUDOTIME_MIN_TOTAL_CELLS
      )
    )
  }


  lineage_condition_counts <- table(
    factor(
      as.character(
        obj_meta[
          lineage_cells,
          CONDITION_COL,
          drop = TRUE
        ]
      ),
      levels = c(
        CONTROL_GROUP,
        CASE_GROUP
      )
    )
  )


  lineage_replicate_counts <- vapply(
    c(
      CONTROL_GROUP,
      CASE_GROUP
    ),
    function(group_name) {
      length(
        unique(
          as.character(
            obj_meta[
              lineage_cells[
                as.character(
                  obj_meta[
                    lineage_cells,
                    CONDITION_COL,
                    drop = TRUE
                  ]
                ) == group_name
              ],
              BIOLOGICAL_REPLICATE_COL,
              drop = TRUE
            ]
          )
        )
      )
    },
    FUN.VALUE = integer(1)
  )


  formal_support_available <-
    !PSEUDOTIME_FORCE_EXPLORATORY &&
    length(
      lineage_cells
    ) >=
      PSEUDOTIME_DEFAULT_MIN_TOTAL_CELLS &&
    all(
      lineage_condition_counts >=
        PSEUDOTIME_DEFAULT_MIN_CELLS_PER_CONDITION
    ) &&
    all(
      lineage_replicate_counts >=
        MIN_REPLICATES_PER_GROUP
    )


  PSEUDOTIME_ANALYSIS_GRADE <- if (
    formal_support_available
  ) {
    "replicate_supported"
  } else {
    "exploratory_low_support"
  }


  PSEUDOTIME_SUPPORT_LABEL <- if (
    formal_support_available
  ) {
    "Replicate-supported descriptive trajectory"
  } else {
    paste0(
      "EXPLORATORY LOW SUPPORT: insufficient cells or biological replicates ",
      "for a Control-to-Case conclusion"
    )
  }


  if (
    !formal_support_available &&
    !PSEUDOTIME_ALLOW_EXPLORATORY_LOW_SUPPORT
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Profile lacks the configured formal cell/replicate support / ",
        "该谱系未达到正式细胞数或生物学重复门槛。 Total cells: ",
        length(lineage_cells),
        "; condition cells: ",
        paste(
          names(lineage_condition_counts),
          as.integer(lineage_condition_counts),
          sep = "=",
          collapse = "; "
        ),
        "; condition replicates: ",
        paste(
          names(lineage_replicate_counts),
          as.integer(lineage_replicate_counts),
          sep = "=",
          collapse = "; "
        )
      )
    )
  }


  pseudotime_profile_support <- data.frame(
    priority = PSEUDOTIME_PROFILE_PRIORITY,
    profile_id = PSEUDOTIME_PROFILE_ID,
    display_name = PSEUDOTIME_PROFILE_DISPLAY_NAME,
    analysis_grade = PSEUDOTIME_ANALYSIS_GRADE,
    target_cell_types = paste(
      PSEUDOTIME_TARGET_CELL_TYPES,
      collapse = "; "
    ),
    n_total_cells = length(
      lineage_cells
    ),
    n_control_cells = as.integer(
      lineage_condition_counts[[CONTROL_GROUP]]
    ),
    n_case_cells = as.integer(
      lineage_condition_counts[[CASE_GROUP]]
    ),
    n_control_replicates = as.integer(
      lineage_replicate_counts[[CONTROL_GROUP]]
    ),
    n_case_replicates = as.integer(
      lineage_replicate_counts[[CASE_GROUP]]
    ),
    allow_exploratory_low_support =
      PSEUDOTIME_ALLOW_EXPLORATORY_LOW_SUPPORT,
    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    pseudotime_profile_support,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_profile_support.csv"
    )
  )


  if (
    any(
      lineage_condition_counts <
        PSEUDOTIME_MIN_CELLS_PER_CONDITION
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Insufficient target-lineage cells in one or both conditions / ",
        "至少一个condition中的目标谱系细胞不足。 Counts: ",
        paste(
          names(lineage_condition_counts),
          as.integer(lineage_condition_counts),
          sep = "=",
          collapse = "; "
        )
      )
    )
  }


  lineage_obj <- subset(
    obj,
    cells = lineage_cells
  )


  DefaultAssay(
    lineage_obj
  ) <- "RNA"


  # ==========================================================
  # 4. RECOMPUTE WITHIN-LINEAGE REPRESENTATION
  # ==========================================================

  set.seed(
    RANDOM_SEED
  )


  n_hvg_use <- min(
    as.integer(
      PSEUDOTIME_N_HVG
    ),
    nrow(
      lineage_obj
    ) - 1L
  )


  lineage_obj <- FindVariableFeatures(
    lineage_obj,
    assay = "RNA",
    selection.method = "vst",
    nfeatures = n_hvg_use,
    verbose = FALSE
  )


  lineage_hvg <- VariableFeatures(
    lineage_obj,
    assay = "RNA"
  )


  if (
    length(
      lineage_hvg
    ) < 3L
  ) {
    safe_stop_module(
      MODULE,
      "Too few within-lineage variable genes / 谱系内部高变基因不足。"
    )
  }


  # The subset inherits the global RNA scale.data layer, whose feature set may
  # differ from the new within-lineage HVGs. Remove only this derived layer before
  # recomputing it; raw counts and normalized data remain unchanged.
  # 子集会继承全局RNA scale.data，其feature集合可能与新的谱系内部HVG不同。重新
  # ScaleData前只删除这一可再生layer；raw counts和标准化data均保持不变。
  if (
    "scale.data" %in%
    SeuratObject::Layers(
      lineage_obj[["RNA"]]
    )
  ) {
    lineage_obj[["RNA"]]$scale.data <- NULL
  }


  lineage_obj <- ScaleData(
    lineage_obj,
    assay = "RNA",
    features = lineage_hvg,
    verbose = FALSE
  )


  n_pcs_use <- min(
    as.integer(
      PSEUDOTIME_N_PCS
    ),
    length(
      lineage_hvg
    ) - 1L,
    ncol(
      lineage_obj
    ) - 1L
  )


  if (
    n_pcs_use < 2L
  ) {
    safe_stop_module(
      MODULE,
      "Too few cells or features for trajectory PCA / 无法计算拟时PCA。"
    )
  }


  lineage_obj <- RunPCA(
    lineage_obj,
    assay = "RNA",
    features = lineage_hvg,
    npcs = n_pcs_use,
    reduction.name = "pseudotime_pca",
    reduction.key = "PTPC_",
    verbose = FALSE
  )


  lineage_obj <- FindNeighbors(
    lineage_obj,
    reduction = "pseudotime_pca",
    dims = seq_len(
      n_pcs_use
    ),
    graph.name = c(
      "pseudotime_nn",
      "pseudotime_snn"
    ),
    verbose = FALSE
  )


  lineage_obj <- FindClusters(
    lineage_obj,
    graph.name = "pseudotime_snn",
    resolution = PSEUDOTIME_CLUSTER_RESOLUTION,
    random.seed = RANDOM_SEED,
    verbose = FALSE
  )


  lineage_obj$pseudotime_cluster <- as.character(
    Idents(
      lineage_obj
    )
  )


  lineage_obj <- RunUMAP(
    lineage_obj,
    reduction = "pseudotime_pca",
    dims = seq_len(
      n_pcs_use
    ),
    reduction.name = "pseudotime_umap",
    reduction.key = "PTUMAP_",
    n.neighbors = min(
      UMAP_N_NEIGHBORS,
      ncol(lineage_obj) - 1L
    ),
    min.dist = UMAP_MIN_DIST,
    umap.method = "uwot",
    metric = "cosine",
    seed.use = RANDOM_SEED,
    verbose = FALSE
  )


  # ==========================================================
  # 4B. EXPLORATORY WITHIN-LINEAGE SUBCLUSTER MARKERS
  #     探索性谱系内亚群MARKER
  # ==========================================================
  #
  # The same expression-derived subclusters used by Slingshot are described here
  # with one-versus-all Seurat marker screening. This cell-level screen helps
  # characterize state heterogeneity and supports manual labels, but individual
  # cells are not independent biological replicates. Therefore these P values
  # must not be reported as replicate-level condition effects; formal Control/
  # Case inference remains the role of M13 pseudobulk edgeR.
  # 本节使用与Slingshot相同的表达驱动亚群，通过Seurat one-versus-all marker
  # 筛查描述谱系内部的状态异质性，并为人工命名提供依据。但单个cell不是
  # 独立生物学重复，因此这些P值不得解读为重复层面的condition效应；正式
  # Control/Case推断仍由M13 pseudobulk edgeR完成。
  #
  # Numeric subclusters are intentionally preserved. Automatic names based on a
  # few top genes would conceal uncertainty and could propagate an incorrect state
  # label into pseudotime interpretation. The blank annotation template makes the
  # required human review explicit.
  # 数字亚群会被有意保留。若仅根据少数top gene自动命名，会掩盖不确定性，
  # 并可能将错误状态标签传递到拟时解读中。空白annotation template用于明确
  # 记录必须进行的人工审阅。

  Idents(
    lineage_obj
  ) <- factor(
    as.character(
      lineage_obj$pseudotime_cluster
    ),
    levels = sort(
      unique(
        as.character(
          lineage_obj$pseudotime_cluster
        )
      )
    )
  )


  lineage_subcluster_counts <- as.data.frame(
    table(
      pseudotime_cluster = as.character(
        lineage_obj$pseudotime_cluster
      )
    ),
    stringsAsFactors = FALSE
  ) |>
    dplyr::rename(
      n_cells = .data$Freq
    ) |>
    dplyr::mutate(
      pseudotime_cluster = as.character(
        .data$pseudotime_cluster
      )
    )


  lineage_marker_result <- tryCatch(
    Seurat::FindAllMarkers(
      lineage_obj,
      assay = "RNA",
      slot = "data",
      only.pos = TRUE,
      min.pct = LINEAGE_SUBCLUSTER_MARKER_MIN_PCT,
      logfc.threshold = LINEAGE_SUBCLUSTER_MARKER_LOGFC,
      test.use = "wilcox",
      verbose = FALSE
    ),
    error = function(e) {
      warning(
        paste0(
          "Within-lineage marker screening failed / ",
          "谱系内亚群marker筛查失败: ",
          conditionMessage(e)
        ),
        call. = FALSE
      )
      data.frame()
    }
  )


  if (
    nrow(
      lineage_marker_result
    ) > 0L
  ) {
    lineage_marker_result$cluster <- as.character(
      lineage_marker_result$cluster
    )


    marker_logfc_column <- intersect(
      c(
        "avg_log2FC",
        "avg_logFC"
      ),
      colnames(
        lineage_marker_result
      )
    )[[1]]


    lineage_top_markers <- lineage_marker_result |>
      dplyr::group_by(
        .data$cluster
      ) |>
      dplyr::arrange(
        .data$p_val_adj,
        dplyr::desc(
          .data[[marker_logfc_column]]
        ),
        .by_group = TRUE
      ) |>
      dplyr::slice_head(
        n = as.integer(
          LINEAGE_SUBCLUSTER_TOP_MARKERS
        )
      ) |>
      dplyr::ungroup()


    data.table::fwrite(
      lineage_marker_result,
      file.path(
        PSEUDOTIME_TABLE_DIR,
        "lineage_subcluster_markers_all.csv"
      )
    )


    data.table::fwrite(
      lineage_top_markers,
      file.path(
        PSEUDOTIME_TABLE_DIR,
        "lineage_subcluster_top_markers.csv"
      )
    )


    lineage_annotation_template <- lineage_subcluster_counts |>
      dplyr::left_join(
        lineage_top_markers |>
          dplyr::group_by(
            .data$cluster
          ) |>
          dplyr::summarise(
            top_markers = paste(
              .data$gene,
              collapse = "; "
            ),
            .groups = "drop"
          ) |>
          dplyr::rename(
            pseudotime_cluster = .data$cluster
          ),
        by = "pseudotime_cluster"
      ) |>
      dplyr::mutate(
        manual_subcluster_label = "",
        annotation_evidence = "",
        reviewer_note = ""
      )


    data.table::fwrite(
      lineage_annotation_template,
      file.path(
        PSEUDOTIME_TABLE_DIR,
        "lineage_subcluster_annotation_template.csv"
      )
    )


    lineage_marker_panel <- unique(
      lineage_top_markers$gene
    )


    if (
      length(
        lineage_marker_panel
      ) > 0L
    ) {
      p_lineage_marker_dotplot <- Seurat::DotPlot(
        lineage_obj,
        features = lineage_marker_panel,
        assay = "RNA",
        group.by = "pseudotime_cluster",
        scale = TRUE,
        dot.scale = 6,
        cols = c(
          "#2166AC",
          "#B2182B"
        )
      ) +
        ggplot2::labs(
          title = paste0(
            PSEUDOTIME_PROFILE_DISPLAY_NAME,
            ": within-lineage subcluster markers"
          ),
          subtitle = paste0(
            "Top ",
            LINEAGE_SUBCLUSTER_TOP_MARKERS,
            " one-versus-all markers per numeric subcluster; exploratory"
          ),
          x = NULL,
          y = "Subcluster"
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1,
            vjust = 1
          ),
          panel.grid = ggplot2::element_blank(),
          plot.title = ggplot2::element_text(
            face = "bold"
          ),
          legend.position = "right"
        )


      save_plot_transparent(
        file.path(
          PSEUDOTIME_FIGURE_DIR,
          "lineage_subcluster_marker_dotplot.png"
        ),
        p_lineage_marker_dotplot,
        width = max(
          10,
          0.38 * length(
            lineage_marker_panel
          ) + 4
        ),
        height = max(
          6,
          0.55 * nrow(
            lineage_subcluster_counts
          ) + 3
        ),
        dpi = 300
      )


      # Average normalized expression is shown per subcluster rather than one
      # column per cell. This retains all marker classes while avoiding an
      # unreadably wide heatmap dominated by abundant populations.
      # 热图按亚群展示标准化表达的平均值，而不是每个cell占一列。这样既保留
      # 所有marker分类，又避免产生过宽且被高丰度亚群主导的热图。
      lineage_marker_expression <- SeuratObject::LayerData(
        lineage_obj,
        assay = "RNA",
        layer = "data"
      )[
        lineage_marker_panel,
        ,
        drop = FALSE
      ]


      lineage_cluster_levels <- levels(
        Idents(
          lineage_obj
        )
      )


      lineage_average_expression <- vapply(
        lineage_cluster_levels,
        function(cluster_name) {
          cluster_cells <- colnames(
            lineage_obj
          )[
            as.character(
              lineage_obj$pseudotime_cluster
            ) == cluster_name
          ]


          Matrix::rowMeans(
            lineage_marker_expression[
              ,
              cluster_cells,
              drop = FALSE
            ]
          )
        },
        FUN.VALUE = numeric(
          length(
            lineage_marker_panel
          )
        )
      )


      rownames(
        lineage_average_expression
      ) <- lineage_marker_panel


      lineage_average_z <- t(
        scale(
          t(
            lineage_average_expression
          )
        )
      )


      lineage_average_z[
        !is.finite(
          lineage_average_z
        )
      ] <- 0


      marker_source_cluster <- stats::setNames(
        as.character(
          lineage_top_markers$cluster
        ),
        lineage_top_markers$gene
      )[
        lineage_marker_panel
      ]


      lineage_heatmap_long <- as.data.frame(
        lineage_average_z,
        check.names = FALSE
      ) |>
        tibble::rownames_to_column(
          "gene"
        ) |>
        tidyr::pivot_longer(
          cols = -dplyr::all_of(
            "gene"
          ),
          names_to = "pseudotime_cluster",
          values_to = "row_z_score"
        ) |>
        dplyr::mutate(
          source_cluster = marker_source_cluster[
            as.character(
              .data$gene
            )
          ],
          gene = factor(
            .data$gene,
            levels = rev(
              lineage_marker_panel
            )
          ),
          pseudotime_cluster = factor(
            .data$pseudotime_cluster,
            levels = lineage_cluster_levels
          )
        )


      data.table::fwrite(
        lineage_heatmap_long,
        file.path(
          PSEUDOTIME_TABLE_DIR,
          "lineage_subcluster_marker_heatmap_matrix.csv"
        )
      )


      p_lineage_marker_heatmap <- ggplot2::ggplot(
        lineage_heatmap_long,
        ggplot2::aes(
          x = .data$pseudotime_cluster,
          y = .data$gene,
          fill = .data$row_z_score
        )
      ) +
        ggplot2::geom_tile(
          colour = "white",
          linewidth = 0.25
        ) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(
            source_cluster
          ),
          scales = "free_y",
          space = "free_y",
          switch = "y"
        ) +
        ggplot2::scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          name = "Row Z-score"
        ) +
        ggplot2::labs(
          title = paste0(
            PSEUDOTIME_PROFILE_DISPLAY_NAME,
            ": classified subcluster marker heatmap"
          ),
          subtitle = paste0(
            "Rows split by marker-source subcluster; columns are within-lineage ",
            "numeric subclusters"
          ),
          x = "Within-lineage subcluster",
          y = NULL
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          strip.placement = "outside",
          strip.text.y.left = ggplot2::element_text(
            angle = 0,
            face = "bold"
          ),
          axis.text.y = ggplot2::element_text(
            size = 8
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          )
        )


      save_plot_transparent(
        file.path(
          PSEUDOTIME_FIGURE_DIR,
          "lineage_subcluster_marker_heatmap.png"
        ),
        p_lineage_marker_heatmap,
        width = max(
          8,
          0.70 * length(
            lineage_cluster_levels
          ) + 4
        ),
        height = max(
          8,
          0.25 * length(
            lineage_marker_panel
          ) + 3
        ),
        dpi = 300
      )
    }
  } else {
    data.table::fwrite(
      lineage_subcluster_counts |>
        dplyr::mutate(
          top_markers = NA_character_,
          manual_subcluster_label = "",
          annotation_evidence = "",
          reviewer_note = ""
        ),
      file.path(
        PSEUDOTIME_TABLE_DIR,
        "lineage_subcluster_annotation_template.csv"
      )
    )
  }


  # ==========================================================
  # 5. ROOT AND TERMINAL SUBCLUSTER SELECTION
  # ==========================================================

  trajectory_normalized_data <- SeuratObject::LayerData(
    lineage_obj,
    assay = "RNA",
    layer = "data"
  )


  root_marker_genes_present <- intersect(
    PSEUDOTIME_ROOT_MARKER_GENES,
    rownames(
      trajectory_normalized_data
    )
  )


  endpoint_marker_genes_present <- intersect(
    PSEUDOTIME_ENDPOINT_MARKER_GENES,
    rownames(
      trajectory_normalized_data
    )
  )


  # A raw arithmetic mean can be dominated by one abundant transcript and is
  # sensitive to per-cell library complexity. UCell instead ranks genes within
  # each cell and measures enrichment of the complete signature near the top of
  # that ranking. The score is therefore used only for relative comparison among
  # subclusters from this same lineage; it is not a universal biological scale.
  # 简单算术均值可能被单个高表达转录本主导，也容易受每个细胞文库复杂度影响。
  # UCell先在细胞内部对基因排序，再评价整个签名在排名前端的富集程度。本分数只
  # 用于同一谱系内部亚群之间的相对比较，不应解释为跨数据集通用的生物学绝对值。
  ucell_signature_features <- list()


  if (
    length(
      root_marker_genes_present
    ) > 0L
  ) {
    ucell_signature_features[["root_identity"]] <-
      root_marker_genes_present
  }


  if (
    length(
      endpoint_marker_genes_present
    ) > 0L
  ) {
    ucell_signature_features[["endpoint_identity"]] <-
      endpoint_marker_genes_present
  }


  ucell_max_rank_use <- min(
    as.integer(
      PSEUDOTIME_UCELL_MAX_RANK
    ),
    nrow(
      trajectory_normalized_data
    )
  )


  if (
    length(
      ucell_signature_features
    ) > 0L
  ) {
    ucell_scores <- UCell::ScoreSignatures_UCell(
      matrix = trajectory_normalized_data,
      features = ucell_signature_features,
      maxRank = ucell_max_rank_use,
      name = "_UCell",
      missing_genes = "skip",
      ncores = 1
    )
  } else {
    ucell_scores <- data.frame(
      row.names = colnames(
        lineage_obj
      )
    )
  }


  get_ucell_score <- function(
      signature_name
  ) {
    score_column <- paste0(
      signature_name,
      "_UCell"
    )

    if (
      !score_column %in%
      colnames(
        ucell_scores
      )
    ) {
      stats::setNames(
        rep(
          NA_real_,
          ncol(lineage_obj)
        ),
        colnames(lineage_obj)
      )
    } else {
      stats::setNames(
        as.numeric(
          ucell_scores[
            colnames(lineage_obj),
            score_column,
            drop = TRUE
          ]
        ),
        colnames(lineage_obj)
      )
    }
  }


  # The historical metadata names are retained for backward compatibility with
  # downstream tables, but their values are now UCell scores rather than means.
  # 为保持下游表格兼容，沿用原metadata列名；其数值现已是UCell签名分数而非均值。
  lineage_obj$pseudotime_root_marker_score <-
    get_ucell_score(
      "root_identity"
    )[
      colnames(lineage_obj)
    ]


  lineage_obj$pseudotime_endpoint_marker_score <-
    get_ucell_score(
      "endpoint_identity"
    )[
      colnames(lineage_obj)
    ]


  marker_set_audit <- data.frame(
    marker_role = c(
      "root",
      "endpoint"
    ),
    configured_genes = c(
      paste(
        PSEUDOTIME_ROOT_MARKER_GENES,
        collapse = "; "
      ),
      paste(
        PSEUDOTIME_ENDPOINT_MARKER_GENES,
        collapse = "; "
      )
    ),
    genes_present = c(
      paste(
        root_marker_genes_present,
        collapse = "; "
      ),
      paste(
        endpoint_marker_genes_present,
        collapse = "; "
      )
    ),
    n_genes_present = c(
      length(
        root_marker_genes_present
      ),
      length(
        endpoint_marker_genes_present
      )
    ),
    scoring_method = "UCell rank-based signature score",
    ucell_max_rank = ucell_max_rank_use,
    interpretation = c(
      paste0(
        "Higher values support the configured root identity/state; ",
        "relative within-profile evidence only."
      ),
      paste0(
        "Higher values support the configured endpoint identity/state; ",
        "relative within-profile evidence only."
      )
    ),
    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    marker_set_audit,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_endpoint_marker_audit.csv"
    )
  )


  # ----------------------------------------------------------
  # Independent reference-label evidence / 独立参考标签证据
  # ----------------------------------------------------------
  # SingleR compares the current expression profile with the curated ImmGen
  # mouse immune reference. It is deliberately restricted to profiles that
  # provide explicit expected label patterns. A reference label is one auxiliary
  # vote—not a replacement for project-specific annotation or UCell signatures.
  # SingleR将当前表达谱与经过整理的ImmGen小鼠免疫参考进行比较。本步骤仅对明确
  # 配置了预期参考标签模式的profile启用。参考标签只是一项辅助证据，不能替代
  # 本项目注释或UCell签名。
  reference_clusters <- sort(
    unique(
      as.character(
        lineage_obj$pseudotime_cluster
      )
    )
  )


  reference_label_audit <- data.frame(
    pseudotime_cluster = reference_clusters,
    reference_name = "ImmGenData label.main",
    reference_status = "not_requested_for_this_profile",
    predicted_label = NA_character_,
    pruned_label = NA_character_,
    label_used_for_matching = NA_character_,
    delta_next = NA_real_,
    root_label_match = NA_real_,
    endpoint_label_match = NA_real_,
    note = paste0(
      "ImmGen is used only for compatible mouse immune profiles; ",
      "other profiles retain UCell plus condition evidence."
    ),
    stringsAsFactors = FALSE
  )


  reference_requested <- isTRUE(
    PSEUDOTIME_USE_IMMGEN_REFERENCE
  ) &&
    identical(
      species_key,
      "mouse"
    ) &&
    length(
      PSEUDOTIME_REFERENCE_ROOT_LABEL_PATTERNS
    ) > 0L &&
    length(
      PSEUDOTIME_REFERENCE_ENDPOINT_LABEL_PATTERNS
    ) > 0L


  matches_reference_patterns <- function(
      labels,
      patterns
  ) {
    if (
      length(
        patterns
      ) == 0L
    ) {
      return(
        rep(
          FALSE,
          length(labels)
        )
      )
    }

    matched <- Reduce(
      `|`,
      lapply(
        patterns,
        function(pattern) {
          grepl(
            pattern,
            labels,
            ignore.case = TRUE,
            perl = TRUE
          )
        }
      )
    )

    !is.na(labels) & matched
  }


  if (
    reference_requested
  ) {
    reference_packages_available <- all(
      vapply(
        c(
          "SingleR",
          "celldex"
        ),
        requireNamespace,
        quietly = TRUE,
        FUN.VALUE = logical(1)
      )
    )


    if (
      !reference_packages_available
    ) {
      PSEUDOTIME_IMMGEN_REFERENCE_ERROR <- paste0(
        "SingleR and/or celldex is unavailable. Run M00 with ",
        "INSTALL_PACKAGES = TRUE."
      )
    }


    if (
      reference_packages_available &&
      is.null(
        PSEUDOTIME_IMMGEN_REFERENCE
      ) &&
      is.na(
        PSEUDOTIME_IMMGEN_REFERENCE_ERROR
      )
    ) {
      reference_load_result <- tryCatch(
        celldex::ImmGenData(
          ensembl = FALSE
        ),
        error = function(e) {
          structure(
            list(
              message = conditionMessage(e)
            ),
            class = "pseudotime_reference_error"
          )
        }
      )


      if (
        inherits(
          reference_load_result,
          "pseudotime_reference_error"
        )
      ) {
        PSEUDOTIME_IMMGEN_REFERENCE_ERROR <-
          reference_load_result$message
      } else {
        PSEUDOTIME_IMMGEN_REFERENCE <-
          reference_load_result
      }
    }


    if (
      !is.null(
        PSEUDOTIME_IMMGEN_REFERENCE
      )
    ) {
      reference_prediction_result <- tryCatch(
        SingleR::SingleR(
          test = trajectory_normalized_data,
          ref = PSEUDOTIME_IMMGEN_REFERENCE,
          labels = as.character(
            SummarizedExperiment::colData(
              PSEUDOTIME_IMMGEN_REFERENCE
            )$label.main
          ),
          clusters = as.character(
            lineage_obj$pseudotime_cluster
          ),
          assay.type.ref = "logcounts",
          check.missing.test = FALSE,
          num.threads = 1
        ),
        error = function(e) {
          structure(
            list(
              message = conditionMessage(e)
            ),
            class = "pseudotime_reference_error"
          )
        }
      )


      if (
        inherits(
          reference_prediction_result,
          "pseudotime_reference_error"
        )
      ) {
        reference_label_audit$reference_status <-
          "prediction_failed_fallback_used"
        reference_label_audit$note <-
          reference_prediction_result$message
      } else {
        reference_prediction <- as.data.frame(
          reference_prediction_result
        )
        reference_prediction$pseudotime_cluster <-
          rownames(
            reference_prediction
          )


        predicted_labels <- as.character(
          reference_prediction$labels
        )


        pruned_labels <- if (
          "pruned.labels" %in%
          colnames(
            reference_prediction
          )
        ) {
          as.character(
            reference_prediction$pruned.labels
          )
        } else {
          rep(
            NA_character_,
            nrow(reference_prediction)
          )
        }


        labels_for_matching <- ifelse(
          !is.na(pruned_labels) &
            trimws(pruned_labels) != "",
          pruned_labels,
          predicted_labels
        )


        reference_prediction_audit <- data.frame(
          pseudotime_cluster = as.character(
            reference_prediction$pseudotime_cluster
          ),
          predicted_label = predicted_labels,
          pruned_label = pruned_labels,
          label_used_for_matching = labels_for_matching,
          delta_next = if (
            "delta.next" %in%
            colnames(
              reference_prediction
            )
          ) {
            as.numeric(
              reference_prediction$delta.next
            )
          } else {
            NA_real_
          },
          root_label_match = as.numeric(
            matches_reference_patterns(
              labels_for_matching,
              PSEUDOTIME_REFERENCE_ROOT_LABEL_PATTERNS
            )
          ),
          endpoint_label_match = as.numeric(
            matches_reference_patterns(
              labels_for_matching,
              PSEUDOTIME_REFERENCE_ENDPOINT_LABEL_PATTERNS
            )
          ),
          stringsAsFactors = FALSE
        )


        reference_label_audit <- reference_label_audit |>
          dplyr::select(
            dplyr::all_of(
              c(
                "pseudotime_cluster",
                "reference_name"
              )
            )
          ) |>
          dplyr::left_join(
            reference_prediction_audit,
            by = "pseudotime_cluster"
          ) |>
          dplyr::mutate(
            reference_status = "available",
            note = paste0(
              "SingleR cluster-level prediction; pruned label is preferred ",
              "when available."
            )
          )
      }
    } else {
      reference_label_audit$reference_status <-
        "reference_load_failed_fallback_used"
      reference_label_audit$note <- if (
        is.na(
          PSEUDOTIME_IMMGEN_REFERENCE_ERROR
        )
      ) {
        "Unknown ImmGen reference-loading error."
      } else {
        PSEUDOTIME_IMMGEN_REFERENCE_ERROR
      }
    }
  }


  root_reference_available <- reference_requested &&
    any(
      reference_label_audit$root_label_match == 1,
      na.rm = TRUE
    )


  endpoint_reference_available <- reference_requested &&
    any(
      reference_label_audit$endpoint_label_match == 1,
      na.rm = TRUE
    )


  data.table::fwrite(
    reference_label_audit,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_reference_label_audit.csv"
    )
  )

  lineage_meta <- lineage_obj[[]] |>
    tibble::rownames_to_column(
      "cell"
    ) |>
    dplyr::mutate(
      pseudotime_cluster = as.character(
        .data$pseudotime_cluster
      ),
      .trajectory_condition = as.character(
        .data[[CONDITION_COL]]
      ),
      .trajectory_cell_type = as.character(
        .data$cell_type
      ),
      .root_marker_score = as.numeric(
        .data$pseudotime_root_marker_score
      ),
      .endpoint_marker_score = as.numeric(
        .data$pseudotime_endpoint_marker_score
      )
    )


  n_control_replicates_total <- lineage_meta |>
    dplyr::filter(
      .data$.trajectory_condition ==
        CONTROL_GROUP
    ) |>
    dplyr::summarise(
      n = dplyr::n_distinct(
        .data[[BIOLOGICAL_REPLICATE_COL]]
      )
    ) |>
    dplyr::pull(
      "n"
    )


  n_case_replicates_total <- lineage_meta |>
    dplyr::filter(
      .data$.trajectory_condition ==
        CASE_GROUP
    ) |>
    dplyr::summarise(
      n = dplyr::n_distinct(
        .data[[BIOLOGICAL_REPLICATE_COL]]
      )
    ) |>
    dplyr::pull(
      "n"
    )


  minimum_control_replicates <- ceiling(
    PSEUDOTIME_MIN_REPLICATE_COVERAGE *
      n_control_replicates_total
  )


  minimum_case_replicates <- ceiling(
    PSEUDOTIME_MIN_REPLICATE_COVERAGE *
      n_case_replicates_total
  )


  trajectory_cluster_audit <- lineage_meta |>
    dplyr::group_by(
      .data$pseudotime_cluster
    ) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      n_control = sum(
        .data$.trajectory_condition ==
          CONTROL_GROUP
      ),
      n_case = sum(
        .data$.trajectory_condition ==
          CASE_GROUP
      ),
      control_fraction = .data$n_control /
        .data$n_cells,
      case_fraction = .data$n_case /
        .data$n_cells,
      root_identity_fraction = if (
        length(
          PSEUDOTIME_ROOT_CELL_TYPES
        ) > 0L
      ) {
        mean(
          .data$.trajectory_cell_type %in%
            PSEUDOTIME_ROOT_CELL_TYPES
        )
      } else {
        NA_real_
      },
      endpoint_identity_fraction = if (
        length(
          PSEUDOTIME_ENDPOINT_CELL_TYPES
        ) > 0L
      ) {
        mean(
          .data$.trajectory_cell_type %in%
            PSEUDOTIME_ENDPOINT_CELL_TYPES
        )
      } else {
        NA_real_
      },
      root_marker_score = if (
        any(
          is.finite(
            .data$.root_marker_score
          )
        )
      ) {
        mean(
          .data$.root_marker_score,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      endpoint_marker_score = if (
        any(
          is.finite(
            .data$.endpoint_marker_score
          )
        )
      ) {
        mean(
          .data$.endpoint_marker_score,
          na.rm = TRUE
        )
      } else {
        NA_real_
      },
      biological_replicates = dplyr::n_distinct(
        .data[[BIOLOGICAL_REPLICATE_COL]]
      ),
      control_replicates = dplyr::n_distinct(
        .data[[BIOLOGICAL_REPLICATE_COL]][
          .data$.trajectory_condition ==
            CONTROL_GROUP
        ]
      ),
      case_replicates = dplyr::n_distinct(
        .data[[BIOLOGICAL_REPLICATE_COL]][
          .data$.trajectory_condition ==
            CASE_GROUP
        ]
      ),
      .groups = "drop"
    )


  trajectory_cluster_audit <- trajectory_cluster_audit |>
    dplyr::left_join(
      reference_label_audit |>
        dplyr::select(
          dplyr::all_of(
            c(
              "pseudotime_cluster",
              "reference_status",
              "predicted_label",
              "pruned_label",
              "label_used_for_matching",
              "delta_next",
              "root_label_match",
              "endpoint_label_match"
            )
          )
        ),
      by = "pseudotime_cluster"
    )


  signature_is_informative <- function(
      cluster_scores,
      genes_present
  ) {
    finite_scores <- cluster_scores[
      is.finite(
        cluster_scores
      )
    ]

    length(
      genes_present
    ) > 0L &&
      length(
        finite_scores
      ) >= 2L &&
      length(
        unique(
          round(
            finite_scores,
            digits = 12
          )
        )
      ) >= 2L
  }


  root_signature_available <- signature_is_informative(
    trajectory_cluster_audit$root_marker_score,
    root_marker_genes_present
  )


  endpoint_signature_available <- signature_is_informative(
    trajectory_cluster_audit$endpoint_marker_score,
    endpoint_marker_genes_present
  )


  effective_evidence_weights <- function(
      signature_available,
      reference_available
  ) {
    raw_weights <- c(
      signature = if (
        signature_available
      ) {
        PSEUDOTIME_SELECTION_WEIGHT_SIGNATURE
      } else {
        0
      },
      reference = if (
        reference_available
      ) {
        PSEUDOTIME_SELECTION_WEIGHT_REFERENCE
      } else {
        0
      },
      condition = PSEUDOTIME_SELECTION_WEIGHT_CONDITION
    )

    raw_weights /
      sum(
        raw_weights
      )
  }


  root_effective_weights <- effective_evidence_weights(
    root_signature_available,
    root_reference_available
  )


  endpoint_effective_weights <- effective_evidence_weights(
    endpoint_signature_available,
    endpoint_reference_available
  )


  # This table is the compact provenance record for endpoint orientation. It
  # distinguishes configured weights from the effective weights after fallback,
  # making the influence of Control/Case explicit and reproducible.
  # 本表是起终点定向的精简溯源记录，区分预设权重与证据缺失后真正采用的权重，
  # 从而明确、可重复地展示Control/Case在本次profile中发挥了多大作用。
  selection_evidence_audit <- data.frame(
    endpoint_role = c(
      "root",
      "endpoint"
    ),
    signature_method = "UCell rank-based signature score",
    signature_available = c(
      root_signature_available,
      endpoint_signature_available
    ),
    reference_method = "SingleR with ImmGenData label.main",
    reference_requested = reference_requested,
    reference_available = c(
      root_reference_available,
      endpoint_reference_available
    ),
    reference_status = paste(
      unique(
        reference_label_audit$reference_status
      ),
      collapse = "; "
    ),
    condition_anchor = c(
      CONTROL_GROUP,
      CASE_GROUP
    ),
    configured_signature_weight =
      PSEUDOTIME_SELECTION_WEIGHT_SIGNATURE,
    configured_reference_weight =
      PSEUDOTIME_SELECTION_WEIGHT_REFERENCE,
    configured_condition_weight =
      PSEUDOTIME_SELECTION_WEIGHT_CONDITION,
    effective_signature_weight = c(
      root_effective_weights[["signature"]],
      endpoint_effective_weights[["signature"]]
    ),
    effective_reference_weight = c(
      root_effective_weights[["reference"]],
      endpoint_effective_weights[["reference"]]
    ),
    effective_condition_weight = c(
      root_effective_weights[["condition"]],
      endpoint_effective_weights[["condition"]]
    ),
    interpretation = paste0(
      "Weights are renormalized across available evidence; condition remains ",
      "an orientation anchor and is not observed chronological time."
    ),
    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    selection_evidence_audit,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_selection_evidence_audit.csv"
    )
  )


  trajectory_cluster_audit <- trajectory_cluster_audit |>
    dplyr::mutate(
      control_enrichment_rank =
        dplyr::percent_rank(
          .data$control_fraction
        ),
      case_enrichment_rank =
        dplyr::percent_rank(
          .data$case_fraction
        ),
      root_marker_rank = if (
        root_signature_available
      ) {
        dplyr::percent_rank(
          dplyr::coalesce(
            .data$root_marker_score,
            stats::median(
              .data$root_marker_score,
              na.rm = TRUE
            )
          )
        )
      } else {
        0
      },
      endpoint_marker_rank = if (
        endpoint_signature_available
      ) {
        dplyr::percent_rank(
          dplyr::coalesce(
            .data$endpoint_marker_score,
            stats::median(
              .data$endpoint_marker_score,
              na.rm = TRUE
            )
          )
        )
      } else {
        0
      },
      root_reference_score = dplyr::coalesce(
        .data$root_label_match,
        0
      ),
      endpoint_reference_score = dplyr::coalesce(
        .data$endpoint_label_match,
        0
      ),
      root_selection_score =
        root_effective_weights[["signature"]] *
          .data$root_marker_rank +
        root_effective_weights[["reference"]] *
          .data$root_reference_score +
        root_effective_weights[["condition"]] *
          .data$control_enrichment_rank,
      endpoint_selection_score =
        endpoint_effective_weights[["signature"]] *
          .data$endpoint_marker_rank +
        endpoint_effective_weights[["reference"]] *
          .data$endpoint_reference_score +
        endpoint_effective_weights[["condition"]] *
          .data$case_enrichment_rank,
      eligible_cluster = if (
        PSEUDOTIME_ANALYSIS_GRADE ==
          "exploratory_low_support"
      ) {
        .data$n_cells >=
          PSEUDOTIME_MIN_CLUSTER_CELLS
      } else {
        .data$n_cells >=
          PSEUDOTIME_MIN_CLUSTER_CELLS &
          .data$control_replicates >=
            minimum_control_replicates &
          .data$case_replicates >=
            minimum_case_replicates
      }
    )


  eligible_clusters <- trajectory_cluster_audit |>
    dplyr::filter(
      .data$eligible_cluster
    )


  if (
    nrow(
      eligible_clusters
    ) < 2L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Fewer than two sufficiently large within-lineage clusters / ",
        "谱系内部满足细胞数阈值的亚群少于两个。"
      )
    )
  }


  root_candidates <- eligible_clusters


  if (
    length(
      PSEUDOTIME_ROOT_CELL_TYPES
    ) > 0L
  ) {
    root_candidates <- root_candidates |>
      dplyr::filter(
        .data$root_identity_fraction > 0.5
      )
  }


  if (
    nrow(
      root_candidates
    ) == 0L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "No eligible root cluster matched PSEUDOTIME_ROOT_CELL_TYPES / ",
        "没有合格起点亚群满足起点cell type约束。"
      )
    )
  }


  root_candidates <- root_candidates |>
    dplyr::arrange(
      dplyr::desc(
        .data$root_selection_score
      ),
      dplyr::desc(
        .data$control_fraction
      ),
      dplyr::desc(
        dplyr::coalesce(
          .data$root_identity_fraction,
          0
        )
      ),
      dplyr::desc(
        .data$n_cells
      ),
      .data$pseudotime_cluster
    )


  root_cluster <- root_candidates$pseudotime_cluster[[1L]]


  endpoint_candidates <- eligible_clusters |>
    dplyr::filter(
      .data$pseudotime_cluster !=
        root_cluster
    )


  if (
    length(
      PSEUDOTIME_ENDPOINT_CELL_TYPES
    ) > 0L
  ) {
    endpoint_candidates <- endpoint_candidates |>
      dplyr::filter(
        .data$endpoint_identity_fraction > 0.5
      )
  }


  if (
    nrow(
      endpoint_candidates
    ) == 0L
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "No eligible terminal cluster distinct from the root matched the ",
        "configured endpoint / 没有与起点不同且满足终点约束的合格亚群。"
      )
    )
  }


  endpoint_candidates <- endpoint_candidates |>
    dplyr::arrange(
      dplyr::desc(
        .data$endpoint_selection_score
      ),
      dplyr::desc(
        .data$case_fraction
      ),
      dplyr::desc(
        dplyr::coalesce(
          .data$endpoint_identity_fraction,
          0
        )
      ),
      dplyr::desc(
        .data$n_cells
      ),
      .data$pseudotime_cluster
    )


  endpoint_cluster <- endpoint_candidates$pseudotime_cluster[[1L]]


  trajectory_cluster_audit <- trajectory_cluster_audit |>
    dplyr::mutate(
      selected_role = dplyr::case_when(
        .data$pseudotime_cluster ==
          root_cluster ~ "root",
        .data$pseudotime_cluster ==
          endpoint_cluster ~ "endpoint",
        TRUE ~ "intermediate_or_unused"
      )
    )


  data.table::fwrite(
    trajectory_cluster_audit,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_cluster_selection_audit.csv"
    )
  )


  # ==========================================================
  # 6. SLINGSHOT TRAJECTORY IN PCA SPACE
  # ==========================================================

  trajectory_pca <- Embeddings(
    lineage_obj,
    reduction = "pseudotime_pca"
  )[
    ,
    seq_len(
      n_pcs_use
    ),
    drop = FALSE
  ]


  trajectory_cluster_labels <- factor(
    as.character(
      lineage_obj$pseudotime_cluster
    )
  )


  names(
    trajectory_cluster_labels
  ) <- colnames(
    lineage_obj
  )


  slingshot_fit <- slingshot::slingshot(
    trajectory_pca,
    clusterLabels = trajectory_cluster_labels,
    start.clus = root_cluster,
    end.clus = endpoint_cluster,
    stretch = 0
  )


  fitted_lineages <- slingshot::slingLineages(
    slingshot_fit
  )


  if (
    length(
      fitted_lineages
    ) == 0L
  ) {
    safe_stop_module(
      MODULE,
      "Slingshot returned no lineage / Slingshot未返回可用轨迹。"
    )
  }


  lineage_matches <- vapply(
    fitted_lineages,
    function(x) {
      identical(
        as.character(
          x[[1L]]
        ),
        root_cluster
      ) &&
        identical(
          as.character(
            x[[length(x)]]
          ),
          endpoint_cluster
        )
    },
    FUN.VALUE = logical(1)
  )


  selected_lineage <- if (
    any(
      lineage_matches
    )
  ) {
    which(
      lineage_matches
    )[[1L]]
  } else {
    1L
  }


  pseudotime_matrix <- slingshot::slingPseudotime(
    slingshot_fit,
    na = TRUE
  )


  curve_weight_matrix <- slingshot::slingCurveWeights(
    slingshot_fit,
    as.probs = TRUE
  )


  pseudotime_raw <- pseudotime_matrix[
    ,
    selected_lineage
  ]


  curve_weight <- curve_weight_matrix[
    ,
    selected_lineage
  ]


  finite_pseudotime <- is.finite(
    pseudotime_raw
  ) &
    is.finite(
      curve_weight
    ) &
    curve_weight > 0


  if (
    sum(
      finite_pseudotime
    ) <
    PSEUDOTIME_MIN_TOTAL_CELLS
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "Too few cells received finite pseudotime values / ",
        "获得有效拟时值的细胞不足: ",
        sum(finite_pseudotime)
      )
    )
  }


  pseudotime_scaled <- rep(
    NA_real_,
    length(
      pseudotime_raw
    )
  )


  pseudotime_range <- range(
    pseudotime_raw[
      finite_pseudotime
    ]
  )


  if (
    diff(
      pseudotime_range
    ) <= 0
  ) {
    safe_stop_module(
      MODULE,
      "Pseudotime has zero range / 拟时结果没有有效变化范围。"
    )
  }


  pseudotime_scaled[
    finite_pseudotime
  ] <- (
    pseudotime_raw[
      finite_pseudotime
    ] -
      pseudotime_range[[1L]]
  ) /
    diff(
      pseudotime_range
    )


  names(
    pseudotime_scaled
  ) <- rownames(
    trajectory_pca
  )


  names(
    curve_weight
  ) <- rownames(
    trajectory_pca
  )


  lineage_obj$pseudotime <- pseudotime_scaled[
    colnames(
      lineage_obj
    )
  ]


  lineage_obj$pseudotime_curve_weight <- curve_weight[
    colnames(
      lineage_obj
    )
  ]


  lineage_obj$pseudotime_root_cluster <- root_cluster
  lineage_obj$pseudotime_endpoint_cluster <- endpoint_cluster


  saveRDS(
    lineage_obj,
    file.path(
      PSEUDOTIME_OBJECT_DIR,
      "lineage_pseudotime_object.rds"
    )
  )


  # ==========================================================
  # 7. CELL-LEVEL AND REPLICATE-LEVEL TABLES
  # ==========================================================

  umap_coordinates <- Embeddings(
    lineage_obj,
    reduction = "pseudotime_umap"
  )


  pca_coordinates <- Embeddings(
    lineage_obj,
    reduction = "pseudotime_pca"
  )[
    ,
    1:2,
    drop = FALSE
  ]


  trajectory_plot_data <- lineage_obj[[]] |>
    tibble::rownames_to_column(
      "cell"
    ) |>
    dplyr::mutate(
      condition = factor(
        as.character(
          .data[[CONDITION_COL]]
        ),
        levels = c(
          CONTROL_GROUP,
          CASE_GROUP
        )
      ),
      biological_replicate = as.character(
        .data[[BIOLOGICAL_REPLICATE_COL]]
      ),
      cell_type = factor(
        as.character(
          .data$cell_type
        ),
        levels = unique(
          PSEUDOTIME_TARGET_CELL_TYPES
        )
      ),
      pseudotime_cluster = as.character(
        .data$pseudotime_cluster
      ),
      PTUMAP_1 = umap_coordinates[
        .data$cell,
        1
      ],
      PTUMAP_2 = umap_coordinates[
        .data$cell,
        2
      ],
      PTPC_1 = pca_coordinates[
        .data$cell,
        1
      ],
      PTPC_2 = pca_coordinates[
        .data$cell,
        2
      ]
    )


  data.table::fwrite(
    trajectory_plot_data |>
      dplyr::select(
        dplyr::all_of(
          c(
            "cell",
            "biological_replicate",
            "condition",
            "cell_type",
            "pseudotime_cluster",
            "pseudotime",
            "pseudotime_curve_weight",
            "PTUMAP_1",
            "PTUMAP_2",
            "PTPC_1",
            "PTPC_2"
          )
        )
      ),
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_cell_metadata.csv"
    )
  )


  pseudotime_replicate_summary <- trajectory_plot_data |>
    dplyr::filter(
      is.finite(
        .data$pseudotime
      )
    ) |>
    dplyr::group_by(
      .data$biological_replicate,
      .data$condition
    ) |>
    dplyr::summarise(
      n_cells = dplyr::n(),
      mean_pseudotime = mean(
        .data$pseudotime
      ),
      median_pseudotime = stats::median(
        .data$pseudotime
      ),
      q25_pseudotime = stats::quantile(
        .data$pseudotime,
        0.25,
        names = FALSE
      ),
      q75_pseudotime = stats::quantile(
        .data$pseudotime,
        0.75,
        names = FALSE
      ),
      .groups = "drop"
    )


  data.table::fwrite(
    pseudotime_replicate_summary,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_by_biological_replicate.csv"
    )
  )


  condition_medians <- trajectory_plot_data |>
    dplyr::filter(
      is.finite(
        .data$pseudotime
      )
    ) |>
    dplyr::group_by(
      .data$condition
    ) |>
    dplyr::summarise(
      median_pseudotime = stats::median(
        .data$pseudotime
      ),
      .groups = "drop"
    )


  control_median <- condition_medians$median_pseudotime[
    condition_medians$condition ==
      CONTROL_GROUP
  ][[1L]]


  case_median <- condition_medians$median_pseudotime[
    condition_medians$condition ==
      CASE_GROUP
  ][[1L]]


  median_pseudotime_difference <-
    case_median -
    control_median


  direction_consistent <- if (
    PSEUDOTIME_ANALYSIS_GRADE ==
      "replicate_supported"
  ) {
    isTRUE(
      median_pseudotime_difference >=
        PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE
    )
  } else {
    NA
  }


  direction_summary <- data.frame(
    priority = PSEUDOTIME_PROFILE_PRIORITY,
    profile_id = PSEUDOTIME_PROFILE_ID,
    display_name = PSEUDOTIME_PROFILE_DISPLAY_NAME,
    analysis_grade = PSEUDOTIME_ANALYSIS_GRADE,
    control_group = CONTROL_GROUP,
    case_group = CASE_GROUP,
    root_cluster = root_cluster,
    endpoint_cluster = endpoint_cluster,
    root_control_fraction =
      root_candidates$control_fraction[[1L]],
    root_ucell_signature_score =
      root_candidates$root_marker_score[[1L]],
    root_reference_label =
      root_candidates$label_used_for_matching[[1L]],
    root_reference_match =
      root_candidates$root_reference_score[[1L]],
    root_effective_signature_weight =
      root_effective_weights[["signature"]],
    root_effective_reference_weight =
      root_effective_weights[["reference"]],
    root_effective_condition_weight =
      root_effective_weights[["condition"]],
    root_selection_score =
      root_candidates$root_selection_score[[1L]],
    endpoint_case_fraction =
      endpoint_candidates$case_fraction[[1L]],
    endpoint_ucell_signature_score =
      endpoint_candidates$endpoint_marker_score[[1L]],
    endpoint_reference_label =
      endpoint_candidates$label_used_for_matching[[1L]],
    endpoint_reference_match =
      endpoint_candidates$endpoint_reference_score[[1L]],
    endpoint_effective_signature_weight =
      endpoint_effective_weights[["signature"]],
    endpoint_effective_reference_weight =
      endpoint_effective_weights[["reference"]],
    endpoint_effective_condition_weight =
      endpoint_effective_weights[["condition"]],
    endpoint_selection_score =
      endpoint_candidates$endpoint_selection_score[[1L]],
    control_median_pseudotime = control_median,
    case_median_pseudotime = case_median,
    case_minus_control_median_pseudotime = median_pseudotime_difference,
    minimum_direction_difference =
      PSEUDOTIME_DIRECTION_MIN_MEDIAN_DIFFERENCE,
    control_to_case_median_direction_consistent = direction_consistent,
    interpretation = if (
      PSEUDOTIME_ANALYSIS_GRADE ==
        "replicate_supported"
    ) {
      paste0(
        "Root/endpoint orientation combined UCell signatures, compatible ",
        "reference labels and Control/Case enrichment using audited effective ",
        "weights. This is descriptive and not an independent progression test."
      )
    } else {
      paste0(
        "Exploratory low-support trajectory: insufficient cells or biological ",
        "replicates for a Control-to-Case directional conclusion."
      )
    },
    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    direction_summary,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_direction_summary.csv"
    )
  )


  if (
    identical(
      direction_consistent,
      FALSE
    )
  ) {
    warning(
      paste0(
        "The fitted cell-level median pseudotime does not increase from ",
        "CONTROL_GROUP to CASE_GROUP by the configured minimum effect size even ",
        "though root/endpoint clusters were oriented by condition. Interpret the ",
        "trajectory cautiously / 尽管起点和终点按condition定向，但cell层面拟时",
        "中位数增幅未达到设定的最小效应量，请谨慎解释。"
      ),
      call. = FALSE
    )
  }


  if (
    is.na(
      direction_consistent
    )
  ) {
    warning(
      paste0(
        "Exploratory trajectory only: cell or biological-replicate support is ",
        "insufficient for a Control-to-Case direction conclusion / 该轨迹仅用于",
        "探索：细胞数或生物学重复不足，不能判断Control到Case方向。"
      ),
      call. = FALSE
    )
  }


  # ==========================================================
  # 8. TRAJECTORY OVERVIEW FIGURE
  # ==========================================================

  trajectory_condition_colours <- stats::setNames(
    c(
      "#0072B2",
      "#D55E00"
    ),
    c(
      CONTROL_GROUP,
      CASE_GROUP
    )
  )


  target_cell_type_colours <- stats::setNames(
    scales::hue_pal()(
      length(
        unique(
          PSEUDOTIME_TARGET_CELL_TYPES
        )
      )
    ),
    unique(
      PSEUDOTIME_TARGET_CELL_TYPES
    )
  )


  p_trajectory_cell_type <- ggplot2::ggplot(
    trajectory_plot_data,
    ggplot2::aes(
      x = .data$PTUMAP_1,
      y = .data$PTUMAP_2,
      colour = .data$cell_type
    )
  ) +
    ggplot2::geom_point(
      size = 0.55,
      alpha = 0.75
    ) +
    ggplot2::scale_colour_manual(
      values = target_cell_type_colours,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Within-lineage UMAP: cell identity",
      x = "Lineage UMAP 1",
      y = "Lineage UMAP 2",
      colour = "Cell type"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )


  p_trajectory_condition <- ggplot2::ggplot(
    trajectory_plot_data,
    ggplot2::aes(
      x = .data$PTUMAP_1,
      y = .data$PTUMAP_2,
      colour = .data$condition
    )
  ) +
    ggplot2::geom_point(
      size = 0.55,
      alpha = 0.75
    ) +
    ggplot2::scale_colour_manual(
      values = trajectory_condition_colours,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Within-lineage UMAP: condition",
      x = "Lineage UMAP 1",
      y = "Lineage UMAP 2",
      colour = "Condition"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )


  p_trajectory_pseudotime <- ggplot2::ggplot(
    trajectory_plot_data |>
      dplyr::filter(
        is.finite(
          .data$pseudotime
        )
      ),
    ggplot2::aes(
      x = .data$PTUMAP_1,
      y = .data$PTUMAP_2,
      colour = .data$pseudotime
    )
  ) +
    ggplot2::geom_point(
      size = 0.6,
      alpha = 0.85
    ) +
    ggplot2::scale_colour_viridis_c(
      option = "C",
      limits = c(
        0,
        1
      )
    ) +
    ggplot2::labs(
      title = "Inferred relative pseudotime",
      x = "Lineage UMAP 1",
      y = "Lineage UMAP 2",
      colour = "Pseudotime"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )


  selected_curve <- slingshot::slingCurves(
    slingshot_fit
  )[[selected_lineage]]


  curve_plot_data <- data.frame(
    PTPC_1 = selected_curve$s[
      selected_curve$ord,
      1
    ],
    PTPC_2 = selected_curve$s[
      selected_curve$ord,
      2
    ]
  )


  p_trajectory_curve <- ggplot2::ggplot(
    trajectory_plot_data |>
      dplyr::filter(
        is.finite(
          .data$pseudotime
        )
      ),
    ggplot2::aes(
      x = .data$PTPC_1,
      y = .data$PTPC_2,
      colour = .data$pseudotime
    )
  ) +
    ggplot2::geom_point(
      size = 0.5,
      alpha = 0.65
    ) +
    ggplot2::geom_path(
      data = curve_plot_data,
      mapping = ggplot2::aes(
        x = .data$PTPC_1,
        y = .data$PTPC_2
      ),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 1
    ) +
    ggplot2::scale_colour_viridis_c(
      option = "C",
      limits = c(
        0,
        1
      )
    ) +
    ggplot2::labs(
      title = paste0(
        "Slingshot PCA curve: cluster ",
        root_cluster,
        " → ",
        endpoint_cluster
      ),
      x = "Lineage PC 1",
      y = "Lineage PC 2",
      colour = "Pseudotime"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "bottom"
    )


  p_trajectory_overview <- (
    p_trajectory_cell_type |
      p_trajectory_condition
  ) /
    (
      p_trajectory_pseudotime |
        p_trajectory_curve
    ) +
    patchwork::plot_annotation(
      title = paste0(
        "Priority ",
        PSEUDOTIME_PROFILE_PRIORITY,
        ": ",
        PSEUDOTIME_PROFILE_DISPLAY_NAME
      ),
      subtitle = paste0(
        PSEUDOTIME_SUPPORT_LABEL,
        "; ",
        CONTROL_GROUP,
        "-relatively-enriched root → ",
        CASE_GROUP,
        "-relatively-enriched endpoint; inferred state direction, not observed time"
      )
    )


  save_plot_transparent(
    file.path(
      PSEUDOTIME_FIGURE_DIR,
      "pseudotime_trajectory_overview.png"
    ),
    p_trajectory_overview,
    width = 15,
    height = 11,
    dpi = 300
  )


  # ==========================================================
  # 9. CONDITION AND BIOLOGICAL-REPLICATE SUMMARY FIGURE
  # ==========================================================

  p_pseudotime_density <- ggplot2::ggplot(
    trajectory_plot_data |>
      dplyr::filter(
        is.finite(
          .data$pseudotime
        )
      ),
    ggplot2::aes(
      x = .data$pseudotime,
      fill = .data$condition,
      colour = .data$condition
    )
  ) +
    ggplot2::geom_density(
      alpha = 0.22,
      linewidth = 0.8,
      adjust = 1.1
    ) +
    ggplot2::scale_fill_manual(
      values = trajectory_condition_colours,
      drop = FALSE
    ) +
    ggplot2::scale_colour_manual(
      values = trajectory_condition_colours,
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Cell-level pseudotime distribution",
      x = "Relative pseudotime",
      y = "Density",
      fill = "Condition",
      colour = "Condition"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      legend.position = "top"
    )


  p_pseudotime_replicates <- ggplot2::ggplot(
    pseudotime_replicate_summary,
    ggplot2::aes(
      x = .data$biological_replicate,
      y = .data$median_pseudotime,
      colour = .data$condition
    )
  ) +
    ggplot2::geom_errorbar(
      ggplot2::aes(
        ymin = .data$q25_pseudotime,
        ymax = .data$q75_pseudotime
      ),
      width = 0.16,
      linewidth = 0.75
    ) +
    ggplot2::geom_point(
      size = 3.2,
      alpha = 0.95
    ) +
    ggplot2::scale_colour_manual(
      values = trajectory_condition_colours,
      drop = FALSE
    ) +
    ggplot2::facet_grid(
      cols = ggplot2::vars(
        condition
      ),
      scales = "free_x",
      space = "free_x"
    ) +
    ggplot2::coord_cartesian(
      ylim = c(
        0,
        1
      )
    ) +
    ggplot2::labs(
      title = "Biological-replicate median pseudotime",
      x = NULL,
      y = "Median pseudotime (IQR)",
      colour = "Condition"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold"
      ),
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(
        face = "bold"
      ),
      axis.text.x = ggplot2::element_text(
        angle = 35,
        hjust = 1
      ),
      legend.position = "none"
    )


  p_pseudotime_condition_summary <-
    p_pseudotime_density |
    p_pseudotime_replicates +
    patchwork::plot_annotation(
      title = paste0(
        PSEUDOTIME_PROFILE_DISPLAY_NAME,
        ": ",
        CONTROL_GROUP,
        " → ",
        CASE_GROUP
      ),
      subtitle = paste0(
        PSEUDOTIME_SUPPORT_LABEL,
        "; condition helped define direction and is not an independent ",
        "progression test"
      )
    )


  save_plot_transparent(
    file.path(
      PSEUDOTIME_FIGURE_DIR,
      "pseudotime_by_condition_and_replicate.png"
    ),
    p_pseudotime_condition_summary,
    width = 13,
    height = 5.8,
    dpi = 300
  )


  # ==========================================================
  # 10. CONFIGURED FEATURE-GENE TRENDS
  # ==========================================================

  pseudotime_feature_genes_present <- intersect(
    PSEUDOTIME_FEATURE_GENES,
    rownames(
      lineage_obj[["RNA"]]
    )
  )


  pseudotime_feature_genes_missing <- setdiff(
    PSEUDOTIME_FEATURE_GENES,
    pseudotime_feature_genes_present
  )


  if (
    length(
      pseudotime_feature_genes_missing
    ) > 0L
  ) {
    message(
      "Pseudotime feature genes not found / 以下拟时展示基因不存在: ",
      paste(
        pseudotime_feature_genes_missing,
        collapse = ", "
      )
    )
  }


  normalized_lineage_data <- SeuratObject::LayerData(
    lineage_obj,
    assay = "RNA",
    layer = "data"
  )


  valid_pseudotime_cells <- trajectory_plot_data$cell[
    is.finite(
      trajectory_plot_data$pseudotime
    )
  ]


  if (
    length(
      pseudotime_feature_genes_present
    ) > 0L
  ) {

    feature_expression <- as.matrix(
      normalized_lineage_data[
        pseudotime_feature_genes_present,
        valid_pseudotime_cells,
        drop = FALSE
      ]
    )


    feature_trend_data <- as.data.frame(
      t(
        feature_expression
      ),
      check.names = FALSE
    ) |>
      tibble::rownames_to_column(
        "cell"
      ) |>
      tidyr::pivot_longer(
        cols = dplyr::all_of(
          pseudotime_feature_genes_present
        ),
        names_to = "gene",
        values_to = "normalized_expression"
      ) |>
      dplyr::left_join(
        trajectory_plot_data |>
          dplyr::select(
            dplyr::all_of(
              c(
                "cell",
                "pseudotime",
                "condition"
              )
            )
          ),
        by = "cell"
      ) |>
      dplyr::mutate(
        gene = factor(
          .data$gene,
          levels = pseudotime_feature_genes_present
        )
      )


    # Keep every configured gene in the point panels, but fit a GAM only when
    # normalized expression contains at least two distinct finite values. A
    # smooth curve has no information for an all-zero/constant gene and can make
    # mgcv report a numerical fitting warning in sparse low-cell profiles.
    # 所有配置基因都保留在散点分面中；仅当标准化表达至少包含两个不同的有限值时
    # 才拟合GAM。对于全零或恒定表达基因，平滑曲线没有信息量，并可能在细胞较少
    # 的稀疏profile中触发mgcv数值拟合warning。
    feature_trend_smooth_data <- feature_trend_data |>
      dplyr::group_by(
        .data$gene
      ) |>
      dplyr::filter(
        dplyr::n_distinct(
          .data$normalized_expression[
            is.finite(
              .data$normalized_expression
            )
          ]
        ) >= 2L
      ) |>
      dplyr::ungroup()


    p_pseudotime_gene_trends <- ggplot2::ggplot(
      feature_trend_data,
      ggplot2::aes(
        x = .data$pseudotime,
        y = .data$normalized_expression,
        colour = .data$condition
      )
    ) +
      ggplot2::geom_point(
        size = 0.35,
        alpha = 0.18
      ) +
      ggplot2::geom_smooth(
        data = feature_trend_smooth_data,
        mapping = ggplot2::aes(
          group = 1
        ),
        method = "gam",
        formula = y ~ s(
          x,
          bs = "cs",
          k = 5
        ),
        se = FALSE,
        colour = "black",
        linewidth = 0.8
      ) +
      ggplot2::scale_colour_manual(
        values = trajectory_condition_colours,
        drop = FALSE
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(
          gene
        ),
        ncol = 4,
        scales = "free_y"
      ) +
      ggplot2::labs(
        title = paste0(
          PSEUDOTIME_PROFILE_DISPLAY_NAME,
          ": selected genes along inferred pseudotime"
        ),
        subtitle = paste0(
          PSEUDOTIME_SUPPORT_LABEL,
          "; black curve = overall GAM trend"
        ),
        x = "Relative pseudotime",
        y = "Normalized RNA expression",
        colour = "Condition"
      ) +
      ggplot2::theme_classic(base_size = 10) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        strip.text = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "top",
        plot.margin = ggplot2::margin(
          t = 10,
          r = 12,
          b = 10,
          l = 28
        )
      )


    save_plot_transparent(
      file.path(
        PSEUDOTIME_FIGURE_DIR,
        "pseudotime_feature_gene_trends.png"
      ),
      p_pseudotime_gene_trends,
      width = 14,
      height = max(
        5,
        3.1 *
          ceiling(
            length(
              pseudotime_feature_genes_present
            ) / 4
          )
      ),
      dpi = 300
    )
  }


  # ==========================================================
  # 11. EXPLORATORY PSEUDOTIME-ASSOCIATED HEATMAP
  # ==========================================================

  heatmap_candidate_genes <- setdiff(
    intersect(
      lineage_hvg,
      rownames(
        normalized_lineage_data
      )
    ),
    grep(
      paste0(
        "(",
        species_cfg$mito_pattern,
        ")|(",
        species_cfg$ribo_pattern,
        ")"
      ),
      intersect(
        lineage_hvg,
        rownames(
          normalized_lineage_data
        )
      ),
      value = TRUE
    )
  )


  heatmap_expression <- as.matrix(
    normalized_lineage_data[
      heatmap_candidate_genes,
      valid_pseudotime_cells,
      drop = FALSE
    ]
  )


  valid_pseudotime_vector <- trajectory_plot_data$pseudotime[
    match(
      valid_pseudotime_cells,
      trajectory_plot_data$cell
    )
  ]


  pseudotime_correlations <- apply(
    heatmap_expression,
    1,
    function(x) {
      suppressWarnings(
        stats::cor(
          x,
          valid_pseudotime_vector,
          method = "spearman",
          use = "pairwise.complete.obs"
        )
      )
    }
  )


  pseudotime_correlation_table <- data.frame(
    gene = names(
      pseudotime_correlations
    ),
    spearman_rho = as.numeric(
      pseudotime_correlations
    ),
    abs_spearman_rho = abs(
      as.numeric(
        pseudotime_correlations
      )
    ),
    stringsAsFactors = FALSE
  ) |>
    dplyr::filter(
      is.finite(
        .data$spearman_rho
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$abs_spearman_rho
      ),
      .data$gene
    )


  data.table::fwrite(
    pseudotime_correlation_table,
    file.path(
      PSEUDOTIME_TABLE_DIR,
      "pseudotime_dynamic_gene_correlations.csv"
    )
  )


  top_dynamic_genes <- utils::head(
    pseudotime_correlation_table$gene,
    min(
      as.integer(
        PSEUDOTIME_HEATMAP_TOP_GENES
      ),
      nrow(
        pseudotime_correlation_table
      )
    )
  )


  if (
    length(
      top_dynamic_genes
    ) >= 2L
  ) {

    pseudotime_bin <- cut(
      valid_pseudotime_vector,
      breaks = seq(
        0,
        1,
        length.out = as.integer(
          PSEUDOTIME_HEATMAP_BINS
        ) + 1L
      ),
      include.lowest = TRUE,
      labels = FALSE
    )


    observed_bins <- sort(
      unique(
        pseudotime_bin[
          !is.na(
            pseudotime_bin
          )
        ]
      )
    )


    binned_expression <- vapply(
      observed_bins,
      function(bin_id) {
        Matrix::rowMeans(
          normalized_lineage_data[
            top_dynamic_genes,
            valid_pseudotime_cells[
              pseudotime_bin ==
                bin_id
            ],
            drop = FALSE
          ]
        )
      },
      FUN.VALUE = numeric(
        length(
          top_dynamic_genes
        )
      )
    )


    rownames(
      binned_expression
    ) <- top_dynamic_genes


    colnames(
      binned_expression
    ) <- as.character(
      observed_bins
    )


    binned_z <- t(
      scale(
        t(
          binned_expression
        )
      )
    )


    binned_z[
      !is.finite(
        binned_z
      )
    ] <- 0


    gene_peak_bin <- apply(
      binned_z,
      1,
      which.max
    )


    ordered_heatmap_genes <- names(
      sort(
        gene_peak_bin
      )
    )


    heatmap_long <- as.data.frame(
      binned_z[
        ordered_heatmap_genes,
        ,
        drop = FALSE
      ],
      check.names = FALSE
    ) |>
      tibble::rownames_to_column(
        "gene"
      ) |>
      tidyr::pivot_longer(
        cols = -dplyr::all_of(
          "gene"
        ),
        names_to = "bin",
        values_to = "z_score"
      ) |>
      dplyr::mutate(
        gene = factor(
          .data$gene,
          levels = rev(
            ordered_heatmap_genes
          )
        ),
        bin = as.integer(
          .data$bin
        ),
        pseudotime_bin_center = (
          .data$bin -
            0.5
        ) /
          as.integer(
            PSEUDOTIME_HEATMAP_BINS
          )
      )


    data.table::fwrite(
      heatmap_long,
      file.path(
        PSEUDOTIME_TABLE_DIR,
        "pseudotime_heatmap_matrix.csv"
      )
    )


    p_pseudotime_heatmap <- ggplot2::ggplot(
      heatmap_long,
      ggplot2::aes(
        x = .data$pseudotime_bin_center,
        y = .data$gene,
        fill = .data$z_score
      )
    ) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_gradient2(
        low = "#2166AC",
        mid = "white",
        high = "#B2182B",
        midpoint = 0,
        name = "Row Z-score"
      ) +
      ggplot2::labs(
        title = paste0(
          PSEUDOTIME_PROFILE_DISPLAY_NAME,
          ": genes associated with inferred pseudotime"
        ),
        subtitle = paste0(
          PSEUDOTIME_SUPPORT_LABEL,
          "; top ",
          length(
            top_dynamic_genes
          ),
          " HVGs ranked by absolute Spearman correlation; exploratory"
        ),
        x = paste0(
          CONTROL_GROUP,
          "-anchored root  →  ",
          CASE_GROUP,
          "-anchored endpoint"
        ),
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        axis.text.y = ggplot2::element_text(
          size = 8
        ),
        legend.position = "right"
      )


    save_plot_transparent(
      file.path(
        PSEUDOTIME_FIGURE_DIR,
        "pseudotime_dynamic_gene_heatmap.png"
      ),
      p_pseudotime_heatmap,
      width = 10,
      height = max(
        7,
        0.23 *
          length(
            top_dynamic_genes
          ) +
          2.5
      ),
      dpi = 300
    )
  }


  write_module_status(
    MODULE,
    "OK"
  )


  data.frame(
    priority = PSEUDOTIME_PROFILE_PRIORITY,
    profile_id = PSEUDOTIME_PROFILE_ID,
    display_name = PSEUDOTIME_PROFILE_DISPLAY_NAME,
    profile_directory = PSEUDOTIME_PROFILE_DIRECTORY,
    module = BASE_PSEUDOTIME_MODULE,
    status = "OK",
    analysis_grade = PSEUDOTIME_ANALYSIS_GRADE,
    n_cells = length(
      lineage_cells
    ),
    tables_directory = PSEUDOTIME_TABLE_DIR,
    figures_directory = PSEUDOTIME_FIGURE_DIR,
    objects_directory = PSEUDOTIME_OBJECT_DIR,
    message = "",
    stringsAsFactors = FALSE
  )
      },
      error = function(e) {
        error_call <- conditionCall(
          e
        )


        error_message <- paste0(
          conditionMessage(
            e
          ),
          if (
            is.null(
              error_call
            )
          ) {
            ""
          } else {
            paste0(
              " [Call: ",
              paste(
                deparse(
                  error_call
                ),
                collapse = " "
              ),
              "]"
            )
          }
        )


        expected_skip <- grepl(
          paste(
            c(
              "None of the requested profile cell types",
              "Too few target-lineage cells",
              "Insufficient target-lineage cells",
              "lacks the configured formal cell/replicate support",
              "Fewer than two sufficiently large",
              "No eligible root cluster",
              "No eligible terminal cluster"
            ),
            collapse = "|"
          ),
          error_message
        )


        profile_status <- if (
          expected_skip
        ) {
          "SKIPPED"
        } else {
          "FAILED"
        }


        write_module_status(
          MODULE,
          profile_status,
          error_message
        )


        warning(
          paste0(
            "Pseudotime profile ",
            PSEUDOTIME_PROFILE_DISPLAY_NAME,
            " ended with status ",
            profile_status,
            ": ",
            error_message
          ),
          call. = FALSE
        )


        data.frame(
          priority = PSEUDOTIME_PROFILE_PRIORITY,
          profile_id = PSEUDOTIME_PROFILE_ID,
          display_name = PSEUDOTIME_PROFILE_DISPLAY_NAME,
          profile_directory = PSEUDOTIME_PROFILE_DIRECTORY,
          module = BASE_PSEUDOTIME_MODULE,
          status = profile_status,
          analysis_grade = PSEUDOTIME_ANALYSIS_GRADE,
          n_cells = length(
            lineage_cells
          ),
          tables_directory = PSEUDOTIME_TABLE_DIR,
          figures_directory = PSEUDOTIME_FIGURE_DIR,
          objects_directory = PSEUDOTIME_OBJECT_DIR,
          message = error_message,
          stringsAsFactors = FALSE
        )
      }
    )


    pseudotime_profile_status_list[[
      PSEUDOTIME_PROFILE_ID
    ]] <- pseudotime_profile_outcome
  }


  pseudotime_profile_status <- dplyr::bind_rows(
    pseudotime_profile_status_list
  ) |>
    dplyr::arrange(
      .data$priority
    )


  data.table::fwrite(
    pseudotime_profile_status,
    file.path(
      BASE_PSEUDOTIME_MDIR,
      "tables",
      "pseudotime_profile_status.csv"
    )
  )


  successful_profile_tables <- pseudotime_profile_status$tables_directory[
    pseudotime_profile_status$status ==
      "OK"
  ]


  combined_direction_summary <- dplyr::bind_rows(
    lapply(
      successful_profile_tables,
      function(profile_table_directory) {
        direction_file <- file.path(
          profile_table_directory,
          "pseudotime_direction_summary.csv"
        )

        if (
          file.exists(
            direction_file
          )
        ) {
          data.table::fread(
            direction_file
          )
        } else {
          NULL
        }
      }
    )
  )


  if (
    nrow(
      combined_direction_summary
    ) > 0L
  ) {
    data.table::fwrite(
      combined_direction_summary,
      file.path(
        BASE_PSEUDOTIME_MDIR,
        "tables",
        "combined_pseudotime_direction_summary.csv"
      )
    )
  }


  if (
    length(
      successful_profile_tables
    ) == 0L
  ) {
    safe_stop_module(
      BASE_PSEUDOTIME_MODULE,
      paste0(
        "No configured pseudotime profile completed successfully / ",
        "没有预设谱系成功完成。"
      )
    )
  }


  write_module_status(
    BASE_PSEUDOTIME_MODULE,
    "OK",
    paste0(
      length(
        successful_profile_tables
      ),
      "/",
      length(
        PSEUDOTIME_LINEAGE_PROFILES
      ),
      " profiles completed"
    )
  )
}



# ============================================================
# M14C. EXPLORATORY CELL-CELL COMMUNICATION / 探索性细胞互作
# ============================================================
#
# 【Module scope / 模块范围】
# CellChat combines normalized ligand/receptor expression with a curated prior
# interaction database to nominate communication patterns among annotated cell
# types. Separate networks are fitted for the configured CONTROL_GROUP and
# CASE_GROUP. The comparison uses only cell types that meet cell-count and
# biological-replicate support thresholds in both conditions.
# CellChat将标准化配体/受体表达与经整理的先验互作数据库结合，用于提出
# 已注釄cell type之间的候选通信模式。模块对配置中的CONTROL_GROUP和CASE_GROUP
# 分别建网；组间比较只使用在两组中均达到细胞数及生物学重复支持门槛的
# cell type。
#
# 【Replicate support / 重复支持】
# A replicate counts as supporting a cell type only when it contains at least
# CELLCHAT_MIN_CELLS_FOR_INTERACTION cells of that type. CellChat's
# filterCommunication(min.samples=...) then retains candidate interactions found
# in the configured minimum number of independent replicates. This reduces, but
# does not eliminate, sensitivity to unequal cell recovery among samples.
# 只有当某个生物学重复中包含至少CELLCHAT_MIN_CELLS_FOR_INTERACTION个相应
# cell type的细胞时，该重复才计为支持。CellChat的filterCommunication(min.samples=...)
# 进一步要求候选互作在指定最少独立重复中出现。这可降低样本间细胞捕获不均
# 的影响，但不能完全消除。
#
# 【Interpretation boundary / 解读边界】
# Communication probabilities are model scores, not replicate-level effect sizes
# and not evidence of direct physical interaction. A condition difference can
# reflect expression, abundance, composition or technical recovery. Candidate
# pathways should be validated with replicate-aware expression results, spatial
# evidence, protein assays or perturbation experiments where possible.
# 通信概率是模型分数，不是生物学重复层面的效应量，也不是细胞直接物理互作的
# 证据。condition差异可能同时反映表达、丰度、组成或技术捕获差异。候选通路应尽可能
# 结合重复感知的表达结果、空间证据、蛋白实验或干预实验验证。
#
# 【Output / 输出】
#   tables/cellchat_cell_type_support_audit.csv
#   tables/cellchat_interactions_all.csv
#   tables/cellchat_pair_strength_by_condition.csv
#   tables/cellchat_pair_strength_difference.csv
#   tables/cellchat_top_LR_comparison.csv
#   figures/cellchat_network_CONTROL_GROUP.png
#   figures/cellchat_network_CASE_GROUP.png
#   figures/cellchat_interaction_strength_heatmap.png
#   figures/cellchat_interaction_strength_difference.png
#   figures/cellchat_top_LR_comparison.png
#   objects/cellchat_CONTROL_GROUP.rds
#   objects/cellchat_CASE_GROUP.rds
#
# 【Method references / 方法依据】
# Original CellChat framework: https://doi.org/10.1038/s41467-021-21246-9
# Updated protocol: https://doi.org/10.1038/s41596-024-01045-4
# Official maintained implementation: https://github.com/jinworks/CellChat
# ============================================================


if (RUN_M14C_CELLCHAT) {

  MODULE <- "M14C_cell_communication"
  MDIR <- module_dir(MODULE)


  cellchat_object_dir <- file.path(
    MDIR,
    "objects"
  )


  dir.create(
    cellchat_object_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )


  if (
    !requireNamespace(
      "CellChat",
      quietly = TRUE
    )
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "CellChat is not installed / 尚未安装CellChat。\n",
        "Run M00 with INSTALL_PACKAGES = TRUE."
      )
    )
  }


  if (
    is.na(
      species_cfg$cellchat_database
    )
  ) {
    data.table::fwrite(
      data.frame(
        status = "SKIPPED",
        reason = paste0(
          "The installed CellChat database branch is configured only for ",
          "mouse and human; current species = ",
          SPECIES
        ),
        stringsAsFactors = FALSE
      ),
      file.path(
        MDIR,
        "tables",
        "cellchat_module_status.csv"
      )
    )


    write_module_status(
      MODULE,
      "OK",
      "CellChat skipped for a species without a configured database."
    )

  } else {

    obj <- load_checkpoint(
      "M11_annotated"
    )


    validate_complete_annotation(
      obj,
      MODULE
    )


    required_cellchat_meta <- c(
      BIOLOGICAL_REPLICATE_COL,
      CONDITION_COL,
      "cell_type"
    )


    missing_cellchat_meta <- setdiff(
      required_cellchat_meta,
      colnames(
        obj[[]]
      )
    )


    if (
      length(
        missing_cellchat_meta
      ) > 0L
    ) {
      safe_stop_module(
        MODULE,
        paste0(
          "CellChat metadata is missing / CellChat所需metadata缺失: ",
          paste(
            missing_cellchat_meta,
            collapse = ", "
          )
        )
      )
    }


    obj <- join_assay_layers_if_needed(
      obj,
      assay = "RNA"
    )


    if (
      !"RNA" %in%
        SeuratObject::Assays(
          obj
        ) ||
      !"data" %in%
        SeuratObject::Layers(
          obj[["RNA"]]
        )
    ) {
      safe_stop_module(
        MODULE,
        "CellChat requires the normalized RNA data layer / CellChat需要RNA标准化data layer。"
      )
    }


    cellchat_meta <- obj[[]] |>
      tibble::rownames_to_column(
        "cell"
      ) |>
      dplyr::transmute(
        cell = .data$cell,
        cell_type = as.character(
          .data$cell_type
        ),
        condition = as.character(
          .data[[CONDITION_COL]]
        ),
        biological_replicate = as.character(
          .data[[BIOLOGICAL_REPLICATE_COL]]
        )
      ) |>
      dplyr::filter(
        .data$condition %in%
          c(
            CONTROL_GROUP,
            CASE_GROUP
          ),
        !is.na(
          .data$cell_type
        ),
        trimws(
          .data$cell_type
        ) != "",
        !is.na(
          .data$biological_replicate
        ),
        trimws(
          .data$biological_replicate
        ) != ""
      )


    cellchat_support_by_replicate <- cellchat_meta |>
      dplyr::count(
        .data$cell_type,
        .data$condition,
        .data$biological_replicate,
        name = "n_cells_in_replicate"
      ) |>
      dplyr::mutate(
        replicate_supported =
          .data$n_cells_in_replicate >=
            CELLCHAT_MIN_CELLS_FOR_INTERACTION
      )


    cellchat_support_audit <- cellchat_meta |>
      dplyr::count(
        .data$cell_type,
        .data$condition,
        name = "n_cells"
      ) |>
      dplyr::left_join(
        cellchat_support_by_replicate |>
          dplyr::group_by(
            .data$cell_type,
            .data$condition
          ) |>
          dplyr::summarise(
            n_replicates_with_any_cell = dplyr::n_distinct(
              .data$biological_replicate
            ),
            n_supported_replicates = dplyr::n_distinct(
              .data$biological_replicate[
                .data$replicate_supported
              ]
            ),
            .groups = "drop"
          ),
        by = c(
          "cell_type",
          "condition"
        )
      ) |>
      tidyr::complete(
        cell_type,
        condition = c(
          CONTROL_GROUP,
          CASE_GROUP
        ),
        fill = list(
          n_cells = 0L,
          n_replicates_with_any_cell = 0L,
          n_supported_replicates = 0L
        )
      ) |>
      dplyr::group_by(
        .data$cell_type
      ) |>
      dplyr::mutate(
        eligible_in_both_conditions = all(
          .data$n_cells >=
            CELLCHAT_MIN_CELLS_PER_TYPE_CONDITION &
            .data$n_supported_replicates >=
              CELLCHAT_MIN_REPLICATES_PER_TYPE_CONDITION
        ),
        eligibility_reason = dplyr::case_when(
          .data$eligible_in_both_conditions ~ "eligible",
          any(
            .data$n_cells <
              CELLCHAT_MIN_CELLS_PER_TYPE_CONDITION
          ) ~ "insufficient total cells in at least one condition",
          TRUE ~ paste0(
            "fewer than ",
            CELLCHAT_MIN_REPLICATES_PER_TYPE_CONDITION,
            " replicates with >= ",
            CELLCHAT_MIN_CELLS_FOR_INTERACTION,
            " cells in at least one condition"
          )
        )
      ) |>
      dplyr::ungroup() |>
      dplyr::arrange(
        .data$cell_type,
        match(
          .data$condition,
          c(
            CONTROL_GROUP,
            CASE_GROUP
          )
        )
      )


    data.table::fwrite(
      cellchat_support_audit,
      file.path(
        MDIR,
        "tables",
        "cellchat_cell_type_support_audit.csv"
      )
    )


    data.table::fwrite(
      cellchat_support_by_replicate,
      file.path(
        MDIR,
        "tables",
        "cellchat_cell_type_support_by_replicate.csv"
      )
    )


    eligible_cell_types <- cellchat_support_audit |>
      dplyr::filter(
        .data$eligible_in_both_conditions
      ) |>
      dplyr::pull(
        .data$cell_type
      ) |>
      unique() |>
      sort()


    if (
      length(
        eligible_cell_types
      ) < 2L
    ) {
      data.table::fwrite(
        data.frame(
          status = "SKIPPED",
          reason = paste0(
            "Fewer than two cell types met support thresholds in both conditions; ",
            "eligible cell types = ",
            length(
              eligible_cell_types
            )
          ),
          stringsAsFactors = FALSE
        ),
        file.path(
          MDIR,
          "tables",
          "cellchat_module_status.csv"
        )
      )


      write_module_status(
        MODULE,
        "OK",
        "CellChat skipped after support audit; fewer than two common cell types."
      )

    } else {

      cellchat_database_name <- paste0(
        "CellChatDB.",
        species_cfg$cellchat_database
      )


      cellchat_database <- getExportedValue(
        "CellChat",
        cellchat_database_name
      )


      cellchat_database <- CellChat::subsetDB(
        cellchat_database,
        search = CELLCHAT_DATABASE_CATEGORIES,
        key = "annotation"
      )


      cellchat_expression <- SeuratObject::LayerData(
        obj,
        assay = "RNA",
        layer = "data"
      )


      cellchat_results <- list()
      cellchat_interactions <- list()
      cellchat_run_audit <- list()


      # CellChat internally uses future-based operations in several versions.
      # A sequential plan is selected for deterministic, memory-conscious module
      # execution; this changes runtime only, not the communication model.
      # CellChat的若干版本会在内部使用future。本模块选择sequential计划，以提高
      # 可重复性并控制内存；该设置只影响运行时间，不改变通信模型。
      future::plan(
        "sequential"
      )


      for (
        condition_name in
          c(
            CONTROL_GROUP,
            CASE_GROUP
          )
      ) {
        condition_meta <- cellchat_meta |>
          dplyr::filter(
            .data$condition ==
              condition_name,
            .data$cell_type %in%
              eligible_cell_types
          ) |>
          dplyr::mutate(
            cell_type = factor(
              .data$cell_type,
              levels = eligible_cell_types
            ),
            samples = factor(
              .data$biological_replicate
            )
          ) |>
          tibble::column_to_rownames(
            "cell"
          )


        condition_expression <- cellchat_expression[
          ,
          rownames(
            condition_meta
          ),
          drop = FALSE
        ]


        condition_result <- tryCatch(
          {
            cellchat_object <- CellChat::createCellChat(
              object = condition_expression,
              meta = condition_meta,
              group.by = "cell_type",
              datatype = "RNA",
              do.sparse = TRUE
            )


            cellchat_object@DB <- cellchat_database


            cellchat_object <- CellChat::subsetData(
              cellchat_object
            )


            cellchat_object <- CellChat::identifyOverExpressedGenes(
              cellchat_object,
              do.fast = FALSE
            )


            cellchat_object <- CellChat::identifyOverExpressedInteractions(
              cellchat_object
            )


            cellchat_object <- CellChat::computeCommunProb(
              cellchat_object,
              type = "triMean",
              raw.use = TRUE,
              population.size = FALSE,
              nboot = 100,
              seed.use = RANDOM_SEED
            )


            cellchat_object <- CellChat::filterCommunication(
              cellchat_object,
              min.cells = CELLCHAT_MIN_CELLS_FOR_INTERACTION,
              min.samples = CELLCHAT_MIN_REPLICATES_PER_TYPE_CONDITION
            )


            cellchat_object <- CellChat::computeCommunProbPathway(
              cellchat_object
            )


            cellchat_object <- CellChat::aggregateNet(
              cellchat_object
            )


            interaction_table <- CellChat::subsetCommunication(
              cellchat_object,
              thresh = 0.05
            )


            interaction_table$condition <- condition_name


            list(
              object = cellchat_object,
              interactions = interaction_table
            )
          },
          error = function(e) {
            structure(
              list(
                message = conditionMessage(e)
              ),
              class = "cellchat_run_error"
            )
          }
        )


        if (
          inherits(
            condition_result,
            "cellchat_run_error"
          )
        ) {
          cellchat_run_audit[[condition_name]] <- data.frame(
            condition = condition_name,
            status = "FAILED",
            n_cells = nrow(
              condition_meta
            ),
            n_interactions = NA_integer_,
            message = condition_result$message,
            stringsAsFactors = FALSE
          )
        } else {
          cellchat_results[[condition_name]] <- condition_result$object
          cellchat_interactions[[condition_name]] <-
            condition_result$interactions


          saveRDS(
            condition_result$object,
            file.path(
              cellchat_object_dir,
              paste0(
                "cellchat_",
                gsub(
                  "[^A-Za-z0-9_-]+",
                  "_",
                  condition_name
                ),
                ".rds"
              )
            )
          )


          cellchat_run_audit[[condition_name]] <- data.frame(
            condition = condition_name,
            status = "OK",
            n_cells = nrow(
              condition_meta
            ),
            n_interactions = nrow(
              condition_result$interactions
            ),
            message = "",
            stringsAsFactors = FALSE
          )
        }
      }


      cellchat_run_audit <- dplyr::bind_rows(
        cellchat_run_audit
      )


      data.table::fwrite(
        cellchat_run_audit,
        file.path(
          MDIR,
          "tables",
          "cellchat_condition_run_audit.csv"
        )
      )


      if (
        any(
          cellchat_run_audit$status !=
            "OK"
        )
      ) {
        safe_stop_module(
          MODULE,
          paste0(
            "CellChat failed in at least one condition / ",
            "CellChat在至少一个condition中失败。 See: ",
            file.path(
              MDIR,
              "tables",
              "cellchat_condition_run_audit.csv"
            )
          )
        )
      }


      interactions_all <- dplyr::bind_rows(
        cellchat_interactions
      ) |>
        dplyr::mutate(
          source = as.character(
            .data$source
          ),
          target = as.character(
            .data$target
          ),
          condition = factor(
            .data$condition,
            levels = c(
              CONTROL_GROUP,
              CASE_GROUP
            )
          )
        )


      data.table::fwrite(
        interactions_all |>
          dplyr::mutate(
            condition = as.character(
              .data$condition
            )
          ),
        file.path(
          MDIR,
          "tables",
          "cellchat_interactions_all.csv"
        )
      )


      pair_strength <- interactions_all |>
        dplyr::group_by(
          .data$condition,
          .data$source,
          .data$target
        ) |>
        dplyr::summarise(
          total_communication_probability = sum(
            .data$prob,
            na.rm = TRUE
          ),
          n_significant_LR_pairs = dplyr::n(),
          .groups = "drop"
        ) |>
        tidyr::complete(
          condition = factor(
            c(
              CONTROL_GROUP,
              CASE_GROUP
            ),
            levels = c(
              CONTROL_GROUP,
              CASE_GROUP
            )
          ),
          source = eligible_cell_types,
          target = eligible_cell_types,
          fill = list(
            total_communication_probability = 0,
            n_significant_LR_pairs = 0L
          )
        )


      data.table::fwrite(
        pair_strength |>
          dplyr::mutate(
            condition = as.character(
              .data$condition
            )
          ),
        file.path(
          MDIR,
          "tables",
          "cellchat_pair_strength_by_condition.csv"
        )
      )


      pair_difference <- pair_strength |>
        dplyr::select(
          "condition",
          "source",
          "target",
          "total_communication_probability",
          "n_significant_LR_pairs"
        ) |>
        tidyr::pivot_wider(
          names_from = "condition",
          values_from = c(
            "total_communication_probability",
            "n_significant_LR_pairs"
          ),
          values_fill = 0
        ) |>
        dplyr::mutate(
          probability_difference =
            .data[[
              paste0(
                "total_communication_probability_",
                CASE_GROUP
              )
            ]] -
              .data[[
                paste0(
                  "total_communication_probability_",
                  CONTROL_GROUP
                )
              ]],
          LR_pair_count_difference =
            .data[[
              paste0(
                "n_significant_LR_pairs_",
                CASE_GROUP
              )
            ]] -
              .data[[
                paste0(
                  "n_significant_LR_pairs_",
                  CONTROL_GROUP
                )
              ]]
        )


      data.table::fwrite(
        pair_difference,
        file.path(
          MDIR,
          "tables",
          "cellchat_pair_strength_difference.csv"
        )
      )


      # ======================================================
      # 1. CIRCULAR COMMUNICATION NETWORKS
      # ======================================================

      cell_type_colours <- stats::setNames(
        scales::hue_pal(
          l = 65,
          c = 100
        )(
          length(
            eligible_cell_types
          )
        ),
        eligible_cell_types
      )


      for (
        condition_name in
          c(
            CONTROL_GROUP,
            CASE_GROUP
          )
      ) {
        condition_edges <- pair_strength |>
          dplyr::filter(
            as.character(
              .data$condition
            ) == condition_name,
            .data$total_communication_probability > 0,
            .data$source !=
              .data$target
          ) |>
          dplyr::transmute(
            from = .data$source,
            to = .data$target,
            weight = .data$total_communication_probability,
            n_LR = .data$n_significant_LR_pairs
          ) |>
          dplyr::slice_max(
            order_by = .data$weight,
            n = as.integer(
              CELLCHAT_NETWORK_TOP_EDGES
            ),
            with_ties = FALSE
          )


        condition_nodes <- cellchat_support_audit |>
          dplyr::filter(
            .data$condition ==
              condition_name,
            .data$cell_type %in%
              eligible_cell_types
          ) |>
          dplyr::transmute(
            name = .data$cell_type,
            n_cells = .data$n_cells
          ) |>
          dplyr::arrange(
            match(
              .data$name,
              eligible_cell_types
            )
          )


        condition_graph <- igraph::graph_from_data_frame(
          condition_edges,
          directed = TRUE,
          vertices = condition_nodes
        )


        p_cellchat_network <- ggraph::ggraph(
          condition_graph,
          layout = "circle"
        ) +
          ggraph::geom_edge_arc(
            ggplot2::aes(
              width = .data$weight,
              alpha = .data$weight
            ),
            colour = "#616161",
            arrow = grid::arrow(
              length = grid::unit(
                2.5,
                "mm"
              ),
              type = "closed"
            ),
            end_cap = ggraph::circle(
              4,
              "mm"
            ),
            strength = 0.18
          ) +
          ggraph::geom_node_point(
            ggplot2::aes(
              size = .data$n_cells,
              colour = .data$name
            )
          ) +
          ggraph::geom_node_text(
            ggplot2::aes(
              label = .data$name
            ),
            repel = TRUE,
            size = 3.2,
            fontface = "bold"
          ) +
          ggraph::scale_edge_width(
            range = c(
              0.25,
              3.5
            ),
            name = "Total\nprobability"
          ) +
          ggraph::scale_edge_alpha(
            range = c(
              0.15,
              0.75
            ),
            guide = "none"
          ) +
          ggplot2::scale_colour_manual(
            values = cell_type_colours,
            guide = "none"
          ) +
          ggplot2::scale_size_continuous(
            range = c(
              5,
              13
            ),
            name = "Cell count"
          ) +
          ggplot2::labs(
            title = paste0(
              "CellChat communication network: ",
              condition_name
            ),
            subtitle = paste0(
              "Top ",
              nrow(
                condition_edges
              ),
              " non-self edges among common supported cell types; edge width ",
              "= summed significant ligand-receptor probability"
            ),
            caption = paste0(
              "Exploratory co-expression model; not direct physical or causal evidence"
            )
          ) +
          ggplot2::theme_void(
            base_size = 10
          ) +
          ggplot2::theme(
            plot.title = ggplot2::element_text(
              face = "bold"
            ),
            legend.position = "right"
          )


        save_plot_transparent(
          file.path(
            MDIR,
            "figures",
            paste0(
              "cellchat_network_",
              gsub(
                "[^A-Za-z0-9_-]+",
                "_",
                condition_name
              ),
              ".png"
            )
          ),
          p_cellchat_network,
          width = 11,
          height = 9,
          dpi = 300
        )
      }


      # ======================================================
      # 2. CONDITION-SPECIFIC AND DIFFERENTIAL HEATMAPS
      # ======================================================

      p_cellchat_strength <- pair_strength |>
        dplyr::mutate(
          source = factor(
            .data$source,
            levels = rev(
              eligible_cell_types
            )
          ),
          target = factor(
            .data$target,
            levels = eligible_cell_types
          )
        ) |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = .data$target,
            y = .data$source,
            fill = log1p(
              .data$total_communication_probability
            )
          )
        ) +
        ggplot2::geom_tile(
          colour = "white",
          linewidth = 0.35
        ) +
        ggplot2::facet_wrap(
          ggplot2::vars(
            condition
          ),
          nrow = 1
        ) +
        ggplot2::scale_fill_gradient(
          low = "white",
          high = "#B2182B",
          name = "log(1 + total\nprobability)"
        ) +
        ggplot2::labs(
          title = "CellChat interaction strength by condition",
          subtitle = "Rows are sender cell types; columns are receiver cell types",
          x = "Receiver",
          y = "Sender"
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1
          ),
          strip.text = ggplot2::element_text(
            face = "bold"
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          )
        )


      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "cellchat_interaction_strength_heatmap.png"
        ),
        p_cellchat_strength,
        width = max(
          12,
          1.0 * length(
            eligible_cell_types
          ) + 6
        ),
        height = max(
          7,
          0.52 * length(
            eligible_cell_types
          ) + 3
        ),
        dpi = 300
      )


      difference_limit <- max(
        abs(
          pair_difference$probability_difference
        ),
        na.rm = TRUE
      )


      if (
        !is.finite(
          difference_limit
        ) ||
        difference_limit == 0
      ) {
        difference_limit <- 1
      }


      p_cellchat_difference <- pair_difference |>
        dplyr::mutate(
          source = factor(
            .data$source,
            levels = rev(
              eligible_cell_types
            )
          ),
          target = factor(
            .data$target,
            levels = eligible_cell_types
          )
        ) |>
        ggplot2::ggplot(
          ggplot2::aes(
            x = .data$target,
            y = .data$source,
            fill = .data$probability_difference
          )
        ) +
        ggplot2::geom_tile(
          colour = "white",
          linewidth = 0.35
        ) +
        ggplot2::scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          limits = c(
            -difference_limit,
            difference_limit
          ),
          name = paste0(
            CASE_GROUP,
            " - ",
            CONTROL_GROUP,
            "\nprobability"
          )
        ) +
        ggplot2::labs(
          title = "CellChat interaction-strength difference",
          subtitle = paste0(
            "Positive = stronger in ",
            CASE_GROUP,
            "; negative = stronger in ",
            CONTROL_GROUP,
            "; descriptive pooled-model contrast"
          ),
          x = "Receiver",
          y = "Sender"
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          )
        )


      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "cellchat_interaction_strength_difference.png"
        ),
        p_cellchat_difference,
        width = max(
          8,
          0.62 * length(
            eligible_cell_types
          ) + 4
        ),
        height = max(
          7,
          0.52 * length(
            eligible_cell_types
          ) + 3
        ),
        dpi = 300
      )


      # ======================================================
      # 3. TOP LIGAND-RECEPTOR COMPARISON
      # ======================================================

      interaction_label_column <- if (
        "interaction_name_2" %in%
          colnames(
            interactions_all
          )
      ) {
        "interaction_name_2"
      } else {
        "interaction_name"
      }


      lr_comparison <- interactions_all |>
        dplyr::transmute(
          condition = as.character(
            .data$condition
          ),
          source = .data$source,
          target = .data$target,
          interaction_name = as.character(
            .data[[interaction_label_column]]
          ),
          pathway_name = as.character(
            .data$pathway_name
          ),
          probability = as.numeric(
            .data$prob
          ),
          p_value = as.numeric(
            .data$pval
          )
        ) |>
        tidyr::complete(
          condition = c(
            CONTROL_GROUP,
            CASE_GROUP
          ),
          nesting(
            source,
            target,
            interaction_name,
            pathway_name
          ),
          fill = list(
            probability = 0,
            p_value = 1
          )
        ) |>
        dplyr::group_by(
          .data$source,
          .data$target,
          .data$interaction_name,
          .data$pathway_name
        ) |>
        dplyr::mutate(
          maximum_probability = max(
            .data$probability,
            na.rm = TRUE
          ),
          probability_difference =
            .data$probability[
              .data$condition ==
                CASE_GROUP
            ] -
              .data$probability[
                .data$condition ==
                  CONTROL_GROUP
              ]
        ) |>
        dplyr::ungroup() |>
        dplyr::group_by(
          .data$source,
          .data$target,
          .data$interaction_name,
          .data$pathway_name
        ) |>
        dplyr::filter(
          dplyr::row_number() == 1L
        ) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          dplyr::desc(
            .data$maximum_probability
          ),
          dplyr::desc(
            abs(
              .data$probability_difference
            )
          )
        ) |>
        dplyr::slice_head(
          n = as.integer(
            CELLCHAT_TOP_LR_PAIRS
          )
        ) |>
        dplyr::select(
          "source",
          "target",
          "interaction_name",
          "pathway_name"
        ) |>
        dplyr::left_join(
          interactions_all |>
            dplyr::transmute(
              condition = as.character(
                .data$condition
              ),
              source = .data$source,
              target = .data$target,
              interaction_name = as.character(
                .data[[interaction_label_column]]
              ),
              pathway_name = as.character(
                .data$pathway_name
              ),
              probability = as.numeric(
                .data$prob
              ),
              p_value = as.numeric(
                .data$pval
              )
            ),
          by = c(
            "source",
            "target",
            "interaction_name",
            "pathway_name"
          )
        ) |>
        tidyr::complete(
          nesting(
            source,
            target,
            interaction_name,
            pathway_name
          ),
          condition = c(
            CONTROL_GROUP,
            CASE_GROUP
          ),
          fill = list(
            probability = 0,
            p_value = 1
          )
        ) |>
        dplyr::mutate(
          comparison_label = paste0(
            .data$source,
            " → ",
            .data$target,
            " | ",
            .data$interaction_name
          ),
          condition = factor(
            .data$condition,
            levels = c(
              CONTROL_GROUP,
              CASE_GROUP
            )
          )
        )


      data.table::fwrite(
        lr_comparison |>
          dplyr::mutate(
            condition = as.character(
              .data$condition
            )
          ),
        file.path(
          MDIR,
          "tables",
          "cellchat_top_LR_comparison.csv"
        )
      )


      lr_label_levels <- lr_comparison |>
        dplyr::group_by(
          .data$comparison_label
        ) |>
        dplyr::summarise(
          maximum_probability = max(
            .data$probability,
            na.rm = TRUE
          ),
          .groups = "drop"
        ) |>
        dplyr::arrange(
          .data$maximum_probability
        ) |>
        dplyr::pull(
          .data$comparison_label
        )


      lr_comparison$comparison_label <- factor(
        lr_comparison$comparison_label,
        levels = lr_label_levels
      )


      p_cellchat_lr <- ggplot2::ggplot(
        lr_comparison,
        ggplot2::aes(
          x = .data$condition,
          y = .data$comparison_label,
          size = .data$probability,
          colour = .data$condition
        )
      ) +
        ggplot2::geom_line(
          ggplot2::aes(
            group = .data$comparison_label
          ),
          colour = "#BDBDBD",
          linewidth = 0.35
        ) +
        ggplot2::geom_point(
          alpha = 0.90
        ) +
        ggplot2::scale_colour_manual(
          values = stats::setNames(
            c(
              "#4DBBD5",
              "#E64B35"
            ),
            c(
              CONTROL_GROUP,
              CASE_GROUP
            )
          ),
          guide = "none"
        ) +
        ggplot2::scale_size_continuous(
          range = c(
            0.5,
            7
          ),
          name = "Communication\nprobability"
        ) +
        ggplot2::labs(
          title = "Top CellChat ligand-receptor comparisons",
          subtitle = paste0(
            "Top ",
            length(
              lr_label_levels
            ),
            " pairs by maximum probability; zero denotes no FDR-significant ",
            "candidate in that condition"
          ),
          x = NULL,
          y = "Sender → receiver | ligand-receptor"
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          panel.grid.major.y = ggplot2::element_line(
            colour = "#EEEEEE",
            linewidth = 0.25
          ),
          panel.grid.minor = ggplot2::element_blank(),
          axis.text.y = ggplot2::element_text(
            size = 7.5
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          )
        )


      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "cellchat_top_LR_comparison.png"
        ),
        p_cellchat_lr,
        width = 12,
        height = max(
          8,
          0.27 * length(
            lr_label_levels
          ) + 3
        ),
        dpi = 300
      )


      data.table::fwrite(
        data.frame(
          status = "OK",
          species_database = cellchat_database_name,
          database_categories = paste(
            CELLCHAT_DATABASE_CATEGORIES,
            collapse = "; "
          ),
          eligible_cell_types = paste(
            eligible_cell_types,
            collapse = "; "
          ),
          n_eligible_cell_types = length(
            eligible_cell_types
          ),
          interpretation = paste0(
            "Exploratory ligand-receptor co-expression model. Condition ",
            "differences are descriptive and do not constitute direct physical ",
            "interaction, causal signalling or replicate-level inference."
          ),
          stringsAsFactors = FALSE
        ),
        file.path(
          MDIR,
          "tables",
          "cellchat_module_status.csv"
        )
      )


      write_module_status(
        MODULE,
        "OK"
      )
    }
  }
}



# ============================================================
# M15. FINAL SAVE / 最终保存
# ============================================================
#
# 【Module scope / 模块范围】
# The final annotated object, evidence synthesis, result manifest and session
# record are assembled for handover and future analysis. Evidence synthesis uses
# deterministic rules to summarize completed modules; it does not refit models,
# modify upstream significance calls or create causal biological conclusions.
# 本模块汇总最终注释对象、证据整合结果、文件清单和会话记录，用于交付和后续分析。
# 证据整合仅按固定规则总结已完成模块，不重新拟合模型、不修改上游显著性判定，也不
# 生成因果性生物学结论。
#
# 【Final acceptance checks / 最终验收检查】
# Confirm that required metadata are complete, raw counts remain accessible, the
# intended PCA/UMAP and identities are present, manual annotation has been written
# back, and the manifest points to all principal tables, figures and checkpoints.
# 应确认必要metadata完整、原始计数仍可访问、预期PCA/UMAP和identity存在、人工注释
# 已写回，并且结果清单包含主要表格、图片和检查点路径。
#
# 【Handover note / 交付说明】
# The final RDS should be accompanied by the user settings, sample metadata,
# package-version log and a concise record of any parameter changes or modules that
# were skipped. The object filename alone is not sufficient provenance.
# 最终RDS应与用户配置、样本metadata、软件版本日志及参数修改或跳过模块的简要记录
# 一并交付；仅提供对象文件名不足以构成完整来源记录。
#
# 【Purpose / 目的】
#
# 保存最终annotated Seurat object，
# 同时记录分析结果文件清单和sessionInfo。
#
#
# 【最终RDS包含什么？】
#
# 至少应包含：
#
# sample_id
# biological_replicate
# condition
# cluster
# cell_type
#
# 同时保留：
#
# raw RNA counts
# normalized expression
# PCA
# UMAP
# cluster identity
#
#
# 【为什么保存最终RDS？】
#
# 后续：
#
# pathway analysis
# gene visualization
# subclustering
# manuscript figures
#
# 都不需要重新运行整个pipeline。
#
# ============================================================


if (RUN_M15_FINAL_SAVE) {

  MODULE <- "M15_final_save"
  MDIR <- module_dir(MODULE)


  obj <- load_checkpoint(
    "M11_annotated"
  )

  validate_complete_annotation(
    obj,
    MODULE
  )


  # ==========================================================
  # 1. FINAL METADATA VALIDATION
  # ==========================================================

  required_final_meta <- c(
    SAMPLE_ID_COL,
    BIOLOGICAL_REPLICATE_COL,
    CONDITION_COL,
    "cluster",
    "cell_type"
  )


  missing_final_meta <- setdiff(
    required_final_meta,
    colnames(
      obj[[]]
    )
  )


  if (
    length(
      missing_final_meta
    ) > 0L
  ) {

    safe_stop_module(
      MODULE,
      paste0(
        "Final object is missing required metadata: ",
        paste(
          missing_final_meta,
          collapse = ", "
        )
      )
    )
  }

  if (!"RNA" %in% SeuratObject::Assays(obj)) {
    safe_stop_module(
      MODULE,
      "Final object is missing the RNA assay / 最终对象缺少RNA assay。"
    )
  }

  rna_layers <- SeuratObject::Layers(
    obj[["RNA"]]
  )

  if (!"data" %in% rna_layers) {
    safe_stop_module(
      MODULE,
      "Final RNA assay is missing the normalized data layer / 最终RNA assay缺少标准化data layer。"
    )
  }

  normalized_data <- tryCatch(
    SeuratObject::LayerData(
      obj,
      assay = "RNA",
      layer = "data"
    ),
    error = function(e) {
      safe_stop_module(
        MODULE,
        paste0(
          "RNA data layer could not be read / RNA data layer无法读取。\n",
          conditionMessage(e)
        )
      )
    }
  )

  expected_dim <- c(
    nrow(obj[["RNA"]]),
    ncol(obj)
  )

  # Compare dimension values rather than storage types. In Seurat v5,
  # dim(LayerData()) may be stored as integer while nrow()/ncol() for an
  # Assay5/Seurat object may be stored as double. identical() would treat
  # c(32285L, 17245L) and c(32285, 17245) as different even though the
  # matrix dimensions are exactly the same.
  #
  # 比较维度数值而不是其底层存储类型。Seurat v5中，dim(LayerData())可能返回
  # integer，而Assay5/Seurat对象的nrow()/ncol()可能返回double。若使用
  # identical()，即使实际行列数完全相同，也会因类型不同而被错误判定为不一致。
  if (
    length(dim(normalized_data)) != 2L ||
    !all(as.numeric(dim(normalized_data)) == expected_dim)
  ) {
    safe_stop_module(
      MODULE,
      paste0(
        "RNA data layer is unreadable or has invalid dimensions / ",
        "RNA data layer不可读或维度错误。"
      )
    )
  }


  # ==========================================================
  # 2. FINAL EVIDENCE SYNTHESIS
  # 2. 最终证据整合
  # ==========================================================
  #
  # This optional subsection reads completed M12-M14C result tables and creates
  # one auditable handover layer. It deliberately separates:
  #
  #   1) replicate-level statistical evidence from M12/M13/M13C;
  #   2) GO/KEGG interpretation derived from significant M13 DEGs;
  #   3) exploratory pseudotime and descriptive CellChat evidence;
  #   4) unavailable or skipped analyses.
  #
  # A pathway name is never treated as proof of a phenotype, and absence of an
  # eligible analysis is never converted into evidence of no biological change.
  # No disease name or cell type is used to choose the primary key cells/genes.
  #
  # 本可选小节读取已完成的M12-M14C结果表，并生成一套可审计的最终交付汇总。
  # 汇总时严格区分：
  #
  #   1）M12/M13/M13C提供的生物学重复层面统计证据；
  #   2）由M13显著DEG派生的GO/KEGG解释；
  #   3）探索性拟时及描述性CellChat证据；
  #   4）不可用或已跳过的分析。
  #
  # 通路名称不会被当作表型证明；未达到分析门槛也不会被错误转换为“无生物学变化”。
  # 关键cell type和关键基因的主排序不使用疾病名称或人工指定的细胞类型。
  # ==========================================================

  if (RUN_FINAL_EVIDENCE_SYNTHESIS) {

    FINAL_REPORT_DIR <- file.path(
      MDIR,
      "report"
    )

    dir.create(
      FINAL_REPORT_DIR,
      recursive = TRUE,
      showWarnings = FALSE
    )


    # --------------------------------------------------------
    # 2.1 SOURCE TABLE AVAILABILITY / 来源表可用性
    # --------------------------------------------------------

    final_source_files <- data.frame(
      source_key = c(
        "composition",
        "pseudobulk_status",
        "significant_DEGs",
        "GO_BP_enrichment",
        "KEGG_enrichment",
        "gene_identifier_mapping",
        "GSVA",
        "pseudotime_profiles",
        "cellchat_support",
        "cellchat_status"
      ),
      evidence_role = c(
        "replicate_level_exploratory_screen",
        "eligibility_audit",
        "replicate_level_statistical_evidence",
        "DEG_dependent_interpretation",
        "DEG_dependent_interpretation",
        "identifier_audit",
        "replicate_level_statistical_evidence",
        "exploratory_trajectory_evidence",
        "descriptive_communication_support",
        "descriptive_communication_status"
      ),
      path = c(
        file.path(
          DIR_MODULES,
          "M12_composition",
          "tables",
          "cell_type_composition_group_comparison_exploratory.csv"
        ),
        file.path(
          DIR_MODULES,
          "M13_pseudobulk_DE",
          "tables",
          "pseudobulk_DE_status.csv"
        ),
        file.path(
          DIR_MODULES,
          "M13_pseudobulk_DE",
          "tables",
          "all_significant_DEGs.csv"
        ),
        file.path(
          DIR_MODULES,
          "M13B_functional_interpretation",
          "tables",
          "GO_BP_enrichment_all.csv"
        ),
        file.path(
          DIR_MODULES,
          "M13B_functional_interpretation",
          "tables",
          "KEGG_enrichment_all.csv"
        ),
        file.path(
          DIR_MODULES,
          "M13B_functional_interpretation",
          "tables",
          "gene_identifier_mapping.csv"
        ),
        file.path(
          DIR_MODULES,
          "M13C_pseudobulk_GSVA",
          "tables",
          "GSVA_all_cell_types_limma.csv"
        ),
        file.path(
          DIR_MODULES,
          "M14B_pseudotime",
          "tables",
          "pseudotime_profile_status.csv"
        ),
        file.path(
          DIR_MODULES,
          "M14C_cell_communication",
          "tables",
          "cellchat_cell_type_support_audit.csv"
        ),
        file.path(
          DIR_MODULES,
          "M14C_cell_communication",
          "tables",
          "cellchat_module_status.csv"
        )
      ),
      stringsAsFactors = FALSE
    )

    final_source_files$available <- file.exists(
      final_source_files$path
    )

    final_source_files$last_modified <- as.character(
      file.info(
        final_source_files$path
      )$mtime
    )

    data.table::fwrite(
      final_source_files,
      file.path(
        MDIR,
        "tables",
        "final_evidence_source_availability.csv"
      )
    )


    final_read_table <- function(source_key) {

      source_row <- final_source_files[
        final_source_files$source_key == source_key,
        ,
        drop = FALSE
      ]

      if (
        nrow(source_row) != 1L ||
        !isTRUE(source_row$available[[1]])
      ) {
        return(NULL)
      }

      tryCatch(
        data.table::fread(
          source_row$path[[1]]
        ),
        error = function(e) {
          warning(
            "Final evidence source could not be read / 最终证据来源表无法读取: ",
            source_row$path[[1]],
            "\n",
            conditionMessage(e),
            call. = FALSE
          )
          NULL
        }
      )
    }


    final_evidence_rules <- data.frame(
      rule_id = c(
        "primary_DEG",
        "primary_GSVA",
        "composition_significant",
        "composition_trend",
        "functional_interpretation",
        "pseudotime",
        "cell_communication",
        "convergent_evidence",
        "supported_evidence",
        "no_signal",
        "insufficient"
      ),
      role = c(
        "statistical",
        "statistical",
        "statistical_exploratory_screen",
        "exploratory",
        "dependent_interpretation",
        "exploratory",
        "descriptive",
        "classification",
        "classification",
        "classification",
        "classification"
      ),
      deterministic_rule = c(
        paste0(
          "M13 pseudobulk PASS; FDR <= ",
          PSEUDOBULK_DE_FDR_THRESHOLD,
          "; abs(logFC) >= ",
          PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ),
        paste0(
          "M13C adj.P.Val <= ",
          GSVA_FDR_THRESHOLD,
          "; abs(pathway logFC) >= ",
          GSVA_MIN_ABS_SCORE_DIFFERENCE
        ),
        paste0(
          "M12 exploratory composition adjusted P < ",
          FINAL_EVIDENCE_FDR_THRESHOLD
        ),
        paste0(
          "M12 nominal P < ",
          FINAL_EVIDENCE_COMPOSITION_TREND_P,
          " but adjusted P >= ",
          FINAL_EVIDENCE_FDR_THRESHOLD
        ),
        paste0(
          "GO/KEGG adjusted P <= ",
          ENRICHMENT_FDR_THRESHOLD,
          "; not counted independently from DEG evidence"
        ),
        "M14B profile retained with its original analysis_grade; exploratory only",
        "M14C eligible common cell type; descriptive pooled CellChat evidence only",
        "At least two significant statistical evidence families; no causal claim",
        "Exactly one significant statistical evidence family",
        "Eligible comparison available but no significant statistical family",
        "No eligible group comparison; absence is not evidence of no change"
      ),
      stringsAsFactors = FALSE
    )

    data.table::fwrite(
      final_evidence_rules,
      file.path(
        MDIR,
        "tables",
        "final_evidence_classification_rules.csv"
      )
    )


    # --------------------------------------------------------
    # 2.2 CELL-TYPE SUPPORT FROM THE FINAL OBJECT
    # 2.2 最终对象中的cell type支持度
    # --------------------------------------------------------

    final_metadata <- obj[[]] |>
      dplyr::transmute(
        cell_type = as.character(
          .data$cell_type
        ),
        condition = as.character(
          .data[[CONDITION_COL]]
        ),
        biological_replicate = as.character(
          .data[[BIOLOGICAL_REPLICATE_COL]]
        )
      )

    final_cell_evidence <- final_metadata |>
      dplyr::group_by(
        .data$cell_type
      ) |>
      dplyr::summarise(
        n_cells_control = sum(
          .data$condition == CONTROL_GROUP,
          na.rm = TRUE
        ),
        n_cells_case = sum(
          .data$condition == CASE_GROUP,
          na.rm = TRUE
        ),
        n_replicates_control = dplyr::n_distinct(
          .data$biological_replicate[
            .data$condition == CONTROL_GROUP
          ]
        ),
        n_replicates_case = dplyr::n_distinct(
          .data$biological_replicate[
            .data$condition == CASE_GROUP
          ]
        ),
        .groups = "drop"
      )


    # --------------------------------------------------------
    # 2.3 PSEUDOBULK ELIGIBILITY AND STATISTICAL KEY GENES
    # 2.3 Pseudobulk可用性及统计学关键基因
    # --------------------------------------------------------

    final_pseudobulk_status <- final_read_table(
      "pseudobulk_status"
    )

    if (!is.null(final_pseudobulk_status)) {
      final_cell_evidence <- final_cell_evidence |>
        dplyr::left_join(
          final_pseudobulk_status |>
            dplyr::transmute(
              cell_type = as.character(
                .data$cell_type
              ),
              pseudobulk_status = as.character(
                .data$status
              ),
              pseudobulk_reason = as.character(
                .data$reason
              ),
              n_control_pseudobulks = as.integer(
                .data$n_control_eligible
              ),
              n_case_pseudobulks = as.integer(
                .data$n_case_eligible
              )
            ),
          by = "cell_type"
        )
    } else {
      final_cell_evidence <- final_cell_evidence |>
        dplyr::mutate(
          pseudobulk_status = "NOT_AVAILABLE",
          pseudobulk_reason = "M13 status table not available",
          n_control_pseudobulks = NA_integer_,
          n_case_pseudobulks = NA_integer_
        )
    }

    final_DEGs <- final_read_table(
      "significant_DEGs"
    )

    final_immediate_early_genes <- c(
      "Fos",
      "Fosb",
      "Fosl1",
      "Fosl2",
      "Jun",
      "Junb",
      "Jund",
      "Egr1",
      "Egr2",
      "Egr3",
      "Dusp1",
      "Ier2",
      "Ier3",
      "Nfkbia",
      "Nfkbiz"
    )

    if (!is.null(final_DEGs) && nrow(final_DEGs) > 0L) {

      final_key_genes <- final_DEGs |>
        dplyr::filter(
          is.finite(
            .data$FDR
          ),
          .data$FDR <=
            PSEUDOBULK_DE_FDR_THRESHOLD,
          abs(
            .data$logFC
          ) >=
            PSEUDOBULK_DE_ABS_LOGFC_THRESHOLD
        ) |>
        dplyr::mutate(
          generic_expression_flag = dplyr::case_when(
            grepl(
              "^mt-",
              .data$gene,
              ignore.case = TRUE
            ) ~ "mitochondrial_gene",
            grepl(
              "^Rp[sl]",
              .data$gene
            ) ~ "ribosomal_gene",
            grepl(
              "^Hsp",
              .data$gene
            ) ~ "heat_shock_gene",
            .data$gene %in%
              final_immediate_early_genes ~
              "immediate_early_or_stress_response_gene",
            TRUE ~ ""
          )
        ) |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::arrange(
          .data$FDR,
          dplyr::desc(
            abs(
              .data$logFC
            )
          ),
          .by_group = TRUE
        ) |>
        dplyr::mutate(
          rank_within_cell_type = dplyr::row_number()
        ) |>
        dplyr::ungroup()

      final_DEG_summary <- final_key_genes |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::summarise(
          n_significant_DEGs = dplyr::n(),
          n_DEGs_higher_in_case = sum(
            .data$logFC > 0,
            na.rm = TRUE
          ),
          n_DEGs_higher_in_control = sum(
            .data$logFC < 0,
            na.rm = TRUE
          ),
          minimum_DEG_FDR = min(
            .data$FDR,
            na.rm = TRUE
          ),
          median_absolute_DEG_logFC = stats::median(
            abs(
              .data$logFC
            ),
            na.rm = TRUE
          ),
          .groups = "drop"
        )

      final_top_key_genes <- final_key_genes |>
        dplyr::filter(
          .data$rank_within_cell_type <=
            as.integer(
              FINAL_EVIDENCE_TOP_GENES_PER_CELL_TYPE
            )
        )

    } else {

      final_key_genes <- data.frame(
        cell_type = character(0),
        contrast = character(0),
        gene = character(0),
        logFC = numeric(0),
        logCPM = numeric(0),
        FDR = numeric(0),
        direction = character(0),
        generic_expression_flag = character(0),
        rank_within_cell_type = integer(0),
        stringsAsFactors = FALSE
      )

      final_top_key_genes <- final_key_genes

      final_DEG_summary <- data.frame(
        cell_type = character(0),
        n_significant_DEGs = integer(0),
        n_DEGs_higher_in_case = integer(0),
        n_DEGs_higher_in_control = integer(0),
        minimum_DEG_FDR = numeric(0),
        median_absolute_DEG_logFC = numeric(0),
        stringsAsFactors = FALSE
      )
    }

    data.table::fwrite(
      final_key_genes,
      file.path(
        MDIR,
        "tables",
        "final_key_genes_statistical.csv"
      )
    )

    data.table::fwrite(
      final_top_key_genes,
      file.path(
        MDIR,
        "tables",
        "final_key_genes_top_for_display.csv"
      )
    )

    final_cell_evidence <- final_cell_evidence |>
      dplyr::left_join(
        final_DEG_summary,
        by = "cell_type"
      )


    # --------------------------------------------------------
    # 2.4 COMPOSITION EVIDENCE / 细胞组成证据
    # --------------------------------------------------------

    final_composition <- final_read_table(
      "composition"
    )

    if (!is.null(final_composition)) {

      final_composition_summary <- final_composition |>
        dplyr::transmute(
          cell_type = as.character(
            .data$cell_type
          ),
          composition_status = as.character(
            .data$status
          ),
          composition_effect_case_minus_control = as.numeric(
            .data$mean_difference_case_minus_control
          ),
          composition_p_value = as.numeric(
            .data$p_value
          ),
          composition_FDR = as.numeric(
            .data$p_adjusted
          ),
          composition_significant = is.finite(
            .data$p_adjusted
          ) &
            .data$p_adjusted <
              FINAL_EVIDENCE_FDR_THRESHOLD,
          composition_exploratory_trend = is.finite(
            .data$p_value
          ) &
            .data$p_value <
              FINAL_EVIDENCE_COMPOSITION_TREND_P &
            (
              !is.finite(
                .data$p_adjusted
              ) |
                .data$p_adjusted >=
                  FINAL_EVIDENCE_FDR_THRESHOLD
            )
        )

      final_cell_evidence <- final_cell_evidence |>
        dplyr::left_join(
          final_composition_summary,
          by = "cell_type"
        )
    } else {
      final_cell_evidence <- final_cell_evidence |>
        dplyr::mutate(
          composition_status = "NOT_AVAILABLE",
          composition_effect_case_minus_control = NA_real_,
          composition_p_value = NA_real_,
          composition_FDR = NA_real_,
          composition_significant = FALSE,
          composition_exploratory_trend = FALSE
        )
    }


    # --------------------------------------------------------
    # 2.5 GO/KEGG DEPENDENT INTERPRETATION
    # 2.5 GO/KEGG依赖性解释
    # --------------------------------------------------------

    final_GO <- final_read_table(
      "GO_BP_enrichment"
    )

    final_KEGG <- final_read_table(
      "KEGG_enrichment"
    )

    final_mapping <- final_read_table(
      "gene_identifier_mapping"
    )

    final_functional_counts <- data.frame(
      cell_type = character(0),
      n_significant_GO_BP_terms = integer(0),
      n_significant_KEGG_pathways = integer(0),
      stringsAsFactors = FALSE
    )

    if (!is.null(final_GO)) {
      final_GO_significant <- final_GO |>
        dplyr::filter(
          is.finite(
            .data$p.adjust
          ),
          .data$p.adjust <=
            ENRICHMENT_FDR_THRESHOLD
        )

      final_GO_counts <- final_GO_significant |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::summarise(
          n_significant_GO_BP_terms = dplyr::n_distinct(
            .data$ID
          ),
          .groups = "drop"
        )
    } else {
      final_GO_significant <- NULL
      final_GO_counts <- data.frame(
        cell_type = character(0),
        n_significant_GO_BP_terms = integer(0)
      )
    }

    if (!is.null(final_KEGG)) {
      final_KEGG_significant <- final_KEGG |>
        dplyr::filter(
          is.finite(
            .data$p.adjust
          ),
          .data$p.adjust <=
            ENRICHMENT_FDR_THRESHOLD
        )

      final_KEGG_counts <- final_KEGG_significant |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::summarise(
          n_significant_KEGG_pathways = dplyr::n_distinct(
            .data$ID
          ),
          .groups = "drop"
        )
    } else {
      final_KEGG_significant <- NULL
      final_KEGG_counts <- data.frame(
        cell_type = character(0),
        n_significant_KEGG_pathways = integer(0)
      )
    }

    final_functional_counts <- dplyr::full_join(
      final_GO_counts,
      final_KEGG_counts,
      by = "cell_type"
    )

    final_cell_evidence <- final_cell_evidence |>
      dplyr::left_join(
        final_functional_counts,
        by = "cell_type"
      )


    final_expand_enrichment_genes <- function(
        enrichment_table,
        ontology_label,
        identifiers_are_entrez = FALSE
    ) {

      empty_result <- data.frame(
        cell_type = character(0),
        direction = character(0),
        ontology = character(0),
        pathway_id = character(0),
        pathway_name = character(0),
        enrichment_FDR = numeric(0),
        gene = character(0),
        stringsAsFactors = FALSE
      )

      if (
        is.null(enrichment_table) ||
        nrow(enrichment_table) == 0L
      ) {
        return(empty_result)
      }

      expanded_rows <- lapply(
        seq_len(
          nrow(enrichment_table)
        ),
        function(i) {

          gene_keys <- strsplit(
            as.character(
              enrichment_table$geneID[[i]]
            ),
            split = "/",
            fixed = TRUE
          )[[1]]

          data.frame(
            cell_type = as.character(
              enrichment_table$cell_type[[i]]
            ),
            direction = as.character(
              enrichment_table$direction[[i]]
            ),
            ontology = ontology_label,
            pathway_id = as.character(
              enrichment_table$ID[[i]]
            ),
            pathway_name = as.character(
              enrichment_table$Description[[i]]
            ),
            enrichment_FDR = as.numeric(
              enrichment_table$p.adjust[[i]]
            ),
            gene_key = gene_keys,
            stringsAsFactors = FALSE
          )
        }
      )

      expanded <- dplyr::bind_rows(
        expanded_rows
      )

      if (identifiers_are_entrez) {

        if (
          is.null(final_mapping) ||
          !all(
            c(
              "SYMBOL",
              "ENTREZID"
            ) %in%
              colnames(
                final_mapping
              )
          )
        ) {
          expanded$gene <- NA_character_
        } else {
          expanded <- expanded |>
            dplyr::mutate(
              gene_key = as.character(
                .data$gene_key
              )
            ) |>
            dplyr::left_join(
              final_mapping |>
                dplyr::transmute(
                  gene_key = as.character(
                    .data$ENTREZID
                  ),
                  gene = as.character(
                    .data$SYMBOL
                  )
                ),
              by = "gene_key"
            )
        }
      } else {
        expanded$gene <- expanded$gene_key
      }

      expanded |>
        dplyr::select(
          "cell_type",
          "direction",
          "ontology",
          "pathway_id",
          "pathway_name",
          "enrichment_FDR",
          "gene"
        ) |>
        dplyr::filter(
          !is.na(
            .data$gene
          ),
          nzchar(
            .data$gene
          )
        )
    }

    final_pathway_driver_genes <- dplyr::bind_rows(
      final_expand_enrichment_genes(
        final_GO_significant,
        "GO Biological Process",
        identifiers_are_entrez = FALSE
      ),
      final_expand_enrichment_genes(
        final_KEGG_significant,
        "KEGG Pathway",
        identifiers_are_entrez = TRUE
      )
    ) |>
      dplyr::distinct() |>
      dplyr::left_join(
        final_key_genes |>
          dplyr::select(
            "cell_type",
            "direction",
            "gene",
            "logFC",
            "FDR",
            "generic_expression_flag",
            "rank_within_cell_type"
          ),
        by = c(
          "cell_type",
          "direction",
          "gene"
        )
      ) |>
      dplyr::arrange(
        .data$cell_type,
        .data$enrichment_FDR,
        .data$FDR,
        dplyr::desc(
          abs(
            .data$logFC
          )
        )
      )

    data.table::fwrite(
      final_pathway_driver_genes,
      file.path(
        MDIR,
        "tables",
        "final_key_genes_pathway_drivers.csv"
      )
    )


    # --------------------------------------------------------
    # 2.6 REPLICATE-AWARE GSVA AND TARGETED KEGG AXES
    # 2.6 生物学重复感知的GSVA及定向KEGG研究轴
    # --------------------------------------------------------

    final_GSVA <- final_read_table(
      "GSVA"
    )

    if (!is.null(final_GSVA)) {

      final_GSVA_summary <- final_GSVA |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::summarise(
          GSVA_tested = TRUE,
          n_significant_GSVA_pathways = sum(
            is.finite(
              .data$adj.P.Val
            ) &
              .data$adj.P.Val <=
                GSVA_FDR_THRESHOLD &
              abs(
                .data$logFC
              ) >=
                GSVA_MIN_ABS_SCORE_DIFFERENCE,
            na.rm = TRUE
          ),
          minimum_GSVA_FDR = if (
            any(
              is.finite(
                .data$adj.P.Val
              )
            )
          ) {
            min(
              .data$adj.P.Val,
              na.rm = TRUE
            )
          } else {
            NA_real_
          },
          maximum_absolute_GSVA_effect = if (
            any(
              is.finite(
                .data$logFC
              )
            )
          ) {
            max(
              abs(
                .data$logFC
              ),
              na.rm = TRUE
            )
          } else {
            NA_real_
          },
          .groups = "drop"
        )

      final_cell_evidence <- final_cell_evidence |>
        dplyr::left_join(
          final_GSVA_summary,
          by = "cell_type"
        )

      final_target_axis_map <- dplyr::bind_rows(
        lapply(
          names(
            FINAL_EVIDENCE_TARGETED_KEGG_AXES
          ),
          function(axis_name) {
            data.frame(
              axis = axis_name,
              pathway_core_id = as.character(
                FINAL_EVIDENCE_TARGETED_KEGG_AXES[[
                  axis_name
                ]]
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      ) |>
        dplyr::mutate(
          axis_label = unname(
            FINAL_EVIDENCE_TARGETED_AXIS_LABELS[
              .data$axis
            ]
          )
        )

      final_targeted_GSVA <- final_GSVA |>
        dplyr::mutate(
          pathway_core_id = sub(
            "^[^0-9]+",
            "",
            as.character(
              .data$pathway_id
            )
          ),
          significant = is.finite(
            .data$adj.P.Val
          ) &
            .data$adj.P.Val <=
              GSVA_FDR_THRESHOLD &
            abs(
              .data$logFC
            ) >=
              GSVA_MIN_ABS_SCORE_DIFFERENCE
        ) |>
        dplyr::inner_join(
          final_target_axis_map,
          by = "pathway_core_id"
        ) |>
        dplyr::mutate(
          evidence_role = "hypothesis_led_secondary_summary",
          interpretation_limit = paste0(
            "Pathway activity association only; does not establish disease, ",
            "metabolite abundance or causal progression."
          )
        ) |>
        dplyr::arrange(
          match(
            .data$axis,
            names(
              FINAL_EVIDENCE_TARGETED_KEGG_AXES
            )
          ),
          .data$adj.P.Val,
          dplyr::desc(
            abs(
              .data$logFC
            )
          )
        )
    } else {

      final_cell_evidence <- final_cell_evidence |>
        dplyr::mutate(
          GSVA_tested = FALSE,
          n_significant_GSVA_pathways = 0L,
          minimum_GSVA_FDR = NA_real_,
          maximum_absolute_GSVA_effect = NA_real_
        )

      final_targeted_GSVA <- data.frame(
        cell_type = character(0),
        pathway_core_id = character(0),
        pathway_id = character(0),
        pathway_name = character(0),
        logFC = numeric(0),
        adj.P.Val = numeric(0),
        significant = logical(0),
        axis = character(0),
        axis_label = character(0),
        evidence_role = character(0),
        interpretation_limit = character(0),
        stringsAsFactors = FALSE
      )
    }

    if (nrow(final_targeted_GSVA) > 0L) {

      final_targeted_GSVA <- final_targeted_GSVA |>
        tidyr::complete(
          cell_type = sort(
            unique(
              final_cell_evidence$cell_type
            )
          ),
          tidyr::nesting(
            axis,
            axis_label,
            pathway_core_id,
            pathway_id,
            pathway_name
          ),
          fill = list(
            significant = FALSE
          )
        ) |>
        dplyr::mutate(
          target_pathway_tested = is.finite(
            .data$logFC
          ),
          evidence_role = ifelse(
            .data$target_pathway_tested,
            "hypothesis_led_secondary_summary",
            "not_tested_or_gene_set_unavailable"
          ),
          interpretation_limit = dplyr::coalesce(
            .data$interpretation_limit,
            paste0(
              "Not tested for this cell type; absence is not evidence of no ",
              "pathway change."
            )
          )
        ) |>
        dplyr::arrange(
          match(
            .data$axis,
            names(
              FINAL_EVIDENCE_TARGETED_KEGG_AXES
            )
          ),
          .data$pathway_name,
          .data$cell_type
        )
    } else {
      final_targeted_GSVA$target_pathway_tested <- logical(0)
    }

    data.table::fwrite(
      final_targeted_GSVA,
      file.path(
        MDIR,
        "tables",
        "final_targeted_KEGG_axis_GSVA_results.csv"
      )
    )


    # --------------------------------------------------------
    # 2.7 EXPLORATORY PSEUDOTIME EVIDENCE
    # 2.7 探索性拟时证据
    # --------------------------------------------------------

    final_pseudotime_support_files <- list.files(
      file.path(
        DIR_MODULES,
        "M14B_pseudotime",
        "tables"
      ),
      pattern = "^pseudotime_profile_support\\.csv$",
      recursive = TRUE,
      full.names = TRUE
    )

    if (length(final_pseudotime_support_files) > 0L) {

      final_pseudotime_support <- data.table::rbindlist(
        lapply(
          final_pseudotime_support_files,
          data.table::fread
        ),
        fill = TRUE
      ) |>
        as.data.frame()

      final_pseudotime_cell_rows <- dplyr::bind_rows(
        lapply(
          seq_len(
            nrow(
              final_pseudotime_support
            )
          ),
          function(i) {

            target_types <- trimws(
              unlist(
                strsplit(
                  as.character(
                    final_pseudotime_support$target_cell_types[[i]]
                  ),
                  split = ";",
                  fixed = TRUE
                )
              )
            )

            data.frame(
              cell_type = target_types,
              profile_id = as.character(
                final_pseudotime_support$profile_id[[i]]
              ),
              display_name = as.character(
                final_pseudotime_support$display_name[[i]]
              ),
              analysis_grade = as.character(
                final_pseudotime_support$analysis_grade[[i]]
              ),
              n_total_cells = as.integer(
                final_pseudotime_support$n_total_cells[[i]]
              ),
              n_control_cells = as.integer(
                final_pseudotime_support$n_control_cells[[i]]
              ),
              n_case_cells = as.integer(
                final_pseudotime_support$n_case_cells[[i]]
              ),
              n_control_replicates = as.integer(
                final_pseudotime_support$n_control_replicates[[i]]
              ),
              n_case_replicates = as.integer(
                final_pseudotime_support$n_case_replicates[[i]]
              ),
              stringsAsFactors = FALSE
            )
          }
        )
      ) |>
        dplyr::filter(
          nzchar(
            .data$cell_type
          )
        )

      final_pseudotime_cell_summary <- final_pseudotime_cell_rows |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::summarise(
          n_pseudotime_profiles = dplyr::n_distinct(
            .data$profile_id
          ),
          pseudotime_analysis_grade = dplyr::case_when(
            any(
              .data$analysis_grade ==
                "replicate_supported"
            ) ~ "replicate_supported_exploratory",
            any(
              .data$analysis_grade ==
                "exploratory_low_support"
            ) ~ "exploratory_low_support",
            TRUE ~ paste(
              sort(
                unique(
                  .data$analysis_grade
                )
              ),
              collapse = "; "
            )
          ),
          pseudotime_profiles = paste(
            sort(
              unique(
                .data$profile_id
              )
            ),
            collapse = "; "
          ),
          .groups = "drop"
        )

      final_cell_evidence <- final_cell_evidence |>
        dplyr::left_join(
          final_pseudotime_cell_summary,
          by = "cell_type"
        )

      data.table::fwrite(
        final_pseudotime_cell_rows,
        file.path(
          MDIR,
          "tables",
          "final_pseudotime_evidence_by_cell_type.csv"
        )
      )
    } else {
      final_cell_evidence <- final_cell_evidence |>
        dplyr::mutate(
          n_pseudotime_profiles = 0L,
          pseudotime_analysis_grade = "NOT_AVAILABLE",
          pseudotime_profiles = ""
        )
    }


    # --------------------------------------------------------
    # 2.8 DESCRIPTIVE CELLCHAT SUPPORT
    # 2.8 描述性CellChat支持
    # --------------------------------------------------------

    final_cellchat_support <- final_read_table(
      "cellchat_support"
    )

    final_cellchat_status <- final_read_table(
      "cellchat_status"
    )

    if (
      !is.null(final_cellchat_support) &&
      nrow(final_cellchat_support) > 0L
    ) {

      final_cellchat_summary <- final_cellchat_support |>
        dplyr::group_by(
          .data$cell_type
        ) |>
        dplyr::summarise(
          communication_eligible = any(
            .data$eligible_in_both_conditions %in%
              TRUE,
            na.rm = TRUE
          ),
          communication_evidence_role = ifelse(
            communication_eligible,
            "descriptive_pooled_model_only",
            "not_eligible"
          ),
          .groups = "drop"
        )

      final_cell_evidence <- final_cell_evidence |>
        dplyr::left_join(
          final_cellchat_summary,
          by = "cell_type"
        )
    } else {
      final_cell_evidence <- final_cell_evidence |>
        dplyr::mutate(
          communication_eligible = FALSE,
          communication_evidence_role = "NOT_AVAILABLE"
        )
    }

    final_cellchat_module_OK <- !is.null(
      final_cellchat_status
    ) &&
      nrow(final_cellchat_status) > 0L &&
      any(
        final_cellchat_status$status == "OK",
        na.rm = TRUE
      )


    # --------------------------------------------------------
    # 2.9 TRANSPARENT EVIDENCE CLASSIFICATION
    # 2.9 透明的证据分类
    # --------------------------------------------------------

    final_cell_evidence <- final_cell_evidence |>
      dplyr::mutate(
        pseudobulk_status = dplyr::coalesce(
          .data$pseudobulk_status,
          "NOT_AVAILABLE"
        ),
        pseudobulk_reason = dplyr::coalesce(
          .data$pseudobulk_reason,
          ""
        ),
        n_significant_DEGs = dplyr::coalesce(
          .data$n_significant_DEGs,
          0L
        ),
        n_DEGs_higher_in_case = dplyr::coalesce(
          .data$n_DEGs_higher_in_case,
          0L
        ),
        n_DEGs_higher_in_control = dplyr::coalesce(
          .data$n_DEGs_higher_in_control,
          0L
        ),
        composition_status = dplyr::coalesce(
          .data$composition_status,
          "NOT_AVAILABLE"
        ),
        composition_significant = dplyr::coalesce(
          .data$composition_significant,
          FALSE
        ),
        composition_exploratory_trend = dplyr::coalesce(
          .data$composition_exploratory_trend,
          FALSE
        ),
        n_significant_GO_BP_terms = dplyr::coalesce(
          .data$n_significant_GO_BP_terms,
          0L
        ),
        n_significant_KEGG_pathways = dplyr::coalesce(
          .data$n_significant_KEGG_pathways,
          0L
        ),
        GSVA_tested = dplyr::coalesce(
          .data$GSVA_tested,
          FALSE
        ),
        n_significant_GSVA_pathways = dplyr::coalesce(
          .data$n_significant_GSVA_pathways,
          0L
        ),
        n_pseudotime_profiles = dplyr::coalesce(
          .data$n_pseudotime_profiles,
          0L
        ),
        pseudotime_analysis_grade = dplyr::coalesce(
          .data$pseudotime_analysis_grade,
          "NOT_AVAILABLE"
        ),
        pseudotime_profiles = dplyr::coalesce(
          .data$pseudotime_profiles,
          ""
        ),
        communication_eligible = dplyr::coalesce(
          .data$communication_eligible,
          FALSE
        ),
        communication_evidence_role = dplyr::coalesce(
          .data$communication_evidence_role,
          "NOT_AVAILABLE"
        ),
        DEG_significant_family =
          .data$pseudobulk_status == "PASS" &
          .data$n_significant_DEGs > 0L,
        GSVA_significant_family =
          .data$GSVA_tested &
          .data$n_significant_GSVA_pathways > 0L,
        composition_significant_family =
          .data$composition_status == "PASS" &
          .data$composition_significant,
        n_significant_statistical_families =
          as.integer(
            .data$DEG_significant_family
          ) +
          as.integer(
            .data$GSVA_significant_family
          ) +
          as.integer(
            .data$composition_significant_family
          ),
        eligible_group_comparison_available =
          .data$pseudobulk_status == "PASS" |
          .data$GSVA_tested |
          .data$composition_status == "PASS",
        exploratory_or_descriptive_evidence =
          .data$composition_exploratory_trend |
          .data$n_pseudotime_profiles > 0L,
        evidence_level = dplyr::case_when(
          .data$n_significant_statistical_families >= 2L ~
            "Convergent statistical evidence",
          .data$n_significant_statistical_families == 1L ~
            "Supported statistical evidence",
          .data$exploratory_or_descriptive_evidence ~
            "Exploratory or descriptive evidence",
          .data$eligible_group_comparison_available ~
            "No significant statistical signal",
          TRUE ~
            "Insufficient or unavailable evidence"
        ),
        evidence_level_order = dplyr::case_when(
          .data$evidence_level ==
            "Convergent statistical evidence" ~ 1L,
          .data$evidence_level ==
            "Supported statistical evidence" ~ 2L,
          .data$evidence_level ==
            "Exploratory or descriptive evidence" ~ 3L,
          .data$evidence_level ==
            "No significant statistical signal" ~ 4L,
          TRUE ~ 5L
        ),
        best_adjusted_P = pmin(
          dplyr::coalesce(
            .data$minimum_DEG_FDR,
            Inf
          ),
          dplyr::coalesce(
            .data$minimum_GSVA_FDR,
            Inf
          ),
          dplyr::coalesce(
            .data$composition_FDR,
            Inf
          )
        ),
        best_adjusted_P = ifelse(
          is.finite(
            .data$best_adjusted_P
          ),
          .data$best_adjusted_P,
          NA_real_
        )
      ) |>
      dplyr::arrange(
        .data$evidence_level_order,
        dplyr::desc(
          .data$n_significant_statistical_families
        ),
        .data$best_adjusted_P,
        dplyr::desc(
          .data$n_significant_DEGs
        ),
        dplyr::desc(
          .data$n_significant_GSVA_pathways
        ),
        .data$cell_type
      ) |>
      dplyr::mutate(
        evidence_priority_rank = dplyr::row_number(),
        rank_basis = paste0(
          "Evidence level; number of significant statistical families; ",
          "best adjusted P; significant DEG count; significant GSVA count. ",
          "No subjective biological weight is used."
        ),
        cellchat_module_status = ifelse(
          final_cellchat_module_OK,
          "OK_descriptive_only",
          "NOT_AVAILABLE_OR_NOT_OK"
        )
      ) |>
      dplyr::select(
        -"evidence_level_order"
      )

    data.table::fwrite(
      final_cell_evidence,
      file.path(
        MDIR,
        "tables",
        "final_key_cell_types_evidence_matrix.csv"
      )
    )


    # --------------------------------------------------------
    # 2.10 FINAL EVIDENCE OVERVIEW FIGURE
    # 2.10 最终证据总览图
    # --------------------------------------------------------

    final_cell_order <- final_cell_evidence$cell_type

    final_evidence_long <- dplyr::bind_rows(
      final_cell_evidence |>
        dplyr::transmute(
          cell_type = .data$cell_type,
          evidence_family = "Composition",
          evidence_status = dplyr::case_when(
            .data$composition_significant ~
              "Statistically significant",
            .data$composition_exploratory_trend ~
              "Exploratory or descriptive",
            .data$composition_status == "PASS" ~
              "No significant signal",
            TRUE ~
              "Unavailable or skipped"
          ),
          display_label = dplyr::case_when(
            .data$composition_significant ~ "FDR",
            .data$composition_exploratory_trend ~ "trend",
            .data$composition_status == "PASS" ~ "0",
            TRUE ~ "NA"
          )
        ),
      final_cell_evidence |>
        dplyr::transmute(
          cell_type = .data$cell_type,
          evidence_family = "Pseudobulk DEG",
          evidence_status = dplyr::case_when(
            .data$DEG_significant_family ~
              "Statistically significant",
            .data$pseudobulk_status == "PASS" ~
              "No significant signal",
            TRUE ~
              "Unavailable or skipped"
          ),
          display_label = dplyr::case_when(
            .data$pseudobulk_status == "PASS" ~
              as.character(
                .data$n_significant_DEGs
              ),
            .data$pseudobulk_status == "SKIPPED" ~ "SKIP",
            TRUE ~ "NA"
          )
        ),
      final_cell_evidence |>
        dplyr::transmute(
          cell_type = .data$cell_type,
          evidence_family = "Pseudobulk GSVA",
          evidence_status = dplyr::case_when(
            .data$GSVA_significant_family ~
              "Statistically significant",
            .data$GSVA_tested ~
              "No significant signal",
            TRUE ~
              "Unavailable or skipped"
          ),
          display_label = ifelse(
            .data$GSVA_tested,
            as.character(
              .data$n_significant_GSVA_pathways
            ),
            "NA"
          )
        ),
      final_cell_evidence |>
        dplyr::transmute(
          cell_type = .data$cell_type,
          evidence_family = "GO / KEGG",
          evidence_status = dplyr::case_when(
            .data$n_significant_GO_BP_terms +
              .data$n_significant_KEGG_pathways > 0L ~
              "DEG-dependent interpretation",
            .data$pseudobulk_status == "PASS" ~
              "No significant signal",
            TRUE ~
              "Unavailable or skipped"
          ),
          display_label = dplyr::case_when(
            .data$n_significant_GO_BP_terms +
              .data$n_significant_KEGG_pathways > 0L ~
              paste0(
                .data$n_significant_GO_BP_terms,
                "/",
                .data$n_significant_KEGG_pathways
              ),
            .data$pseudobulk_status == "PASS" ~ "0/0",
            TRUE ~ "NA"
          )
        ),
      final_cell_evidence |>
        dplyr::transmute(
          cell_type = .data$cell_type,
          evidence_family = "Pseudotime",
          evidence_status = ifelse(
            .data$n_pseudotime_profiles > 0L,
            "Exploratory or descriptive",
            "Unavailable or skipped"
          ),
          display_label = dplyr::case_when(
            grepl(
              "replicate_supported",
              .data$pseudotime_analysis_grade,
              fixed = TRUE
            ) ~ "R",
            .data$n_pseudotime_profiles > 0L ~ "E",
            TRUE ~ "NA"
          )
        ),
      final_cell_evidence |>
        dplyr::transmute(
          cell_type = .data$cell_type,
          evidence_family = "CellChat",
          evidence_status = ifelse(
            .data$communication_eligible,
            "Exploratory or descriptive",
            "Unavailable or skipped"
          ),
          display_label = ifelse(
            .data$communication_eligible,
            "D",
            "NA"
          )
        )
    ) |>
      dplyr::mutate(
        cell_type = factor(
          .data$cell_type,
          levels = rev(
            final_cell_order
          )
        ),
        evidence_family = factor(
          .data$evidence_family,
          levels = c(
            "Composition",
            "Pseudobulk DEG",
            "Pseudobulk GSVA",
            "GO / KEGG",
            "Pseudotime",
            "CellChat"
          )
        )
      )

    data.table::fwrite(
      final_evidence_long |>
        dplyr::mutate(
          cell_type = as.character(
            .data$cell_type
          ),
          evidence_family = as.character(
            .data$evidence_family
          )
        ),
      file.path(
        MDIR,
        "tables",
        "final_cell_type_evidence_plot_data.csv"
      )
    )

    p_final_evidence <- ggplot2::ggplot(
      final_evidence_long,
      ggplot2::aes(
        x = .data$evidence_family,
        y = .data$cell_type,
        fill = .data$evidence_status
      )
    ) +
      ggplot2::geom_tile(
        colour = "white",
        linewidth = 0.45
      ) +
      ggplot2::geom_text(
        ggplot2::aes(
          label = .data$display_label
        ),
        size = 3.1
      ) +
      ggplot2::scale_fill_manual(
        values = c(
          "Statistically significant" = "#B2182B",
          "DEG-dependent interpretation" = "#EF8A62",
          "Exploratory or descriptive" = "#67A9CF",
          "No significant signal" = "#F7F7F7",
          "Unavailable or skipped" = "#BDBDBD"
        ),
        drop = FALSE,
        name = "Evidence status"
      ) +
      ggplot2::labs(
        title = "Final cell-type evidence overview",
        subtitle = paste0(
          "Numbers are significant DEG or pathway counts; GO/KEGG = GO/KEGG ",
          "counts; R = replicate-supported exploratory trajectory; ",
          "E = low-support trajectory; D = descriptive CellChat; ",
          "SKIP/NA are not evidence of no change."
        ),
        x = NULL,
        y = NULL,
        caption = paste0(
          "Primary ranking uses fixed statistical rules only; GO/KEGG, ",
          "pseudotime and CellChat do not add independent statistical evidence."
        )
      ) +
      ggplot2::theme_minimal(
        base_size = 10
      ) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(
          angle = 30,
          hjust = 1
        ),
        plot.title = ggplot2::element_text(
          face = "bold"
        ),
        legend.position = "bottom"
      )

    save_plot_transparent(
      file.path(
        MDIR,
        "figures",
        "final_cell_type_evidence_overview.png"
      ),
      p_final_evidence,
      width = 12,
      height = max(
        7,
        0.40 *
          length(
            final_cell_order
          ) +
          3
      ),
      dpi = 300
    )


    # --------------------------------------------------------
    # 2.11 STATISTICALLY SELECTED KEY-GENE FIGURE
    # 2.11 统计学筛选的关键基因图
    # --------------------------------------------------------

    if (nrow(final_top_key_genes) > 0L) {

      final_top_key_genes_plot <- final_top_key_genes |>
        dplyr::mutate(
          comparison_direction = ifelse(
            .data$logFC > 0,
            paste0(
              "Higher in ",
              CASE_GROUP
            ),
            paste0(
              "Higher in ",
              CONTROL_GROUP
            )
          ),
          negative_log10_FDR = -log10(
            pmax(
              .data$FDR,
              .Machine$double.xmin
            )
          ),
          gene_cell_key = paste(
            .data$gene,
            .data$cell_type,
            sep = "___"
          )
        ) |>
        dplyr::arrange(
          .data$cell_type,
          .data$logFC
        ) |>
        dplyr::mutate(
          gene_cell_key = factor(
            .data$gene_cell_key,
            levels = unique(
              .data$gene_cell_key
            )
          )
        )

      p_final_key_genes <- ggplot2::ggplot(
        final_top_key_genes_plot,
        ggplot2::aes(
          x = .data$logFC,
          y = .data$gene_cell_key,
          colour = .data$comparison_direction
        )
      ) +
        ggplot2::geom_vline(
          xintercept = 0,
          colour = "#BDBDBD",
          linewidth = 0.4
        ) +
        ggplot2::geom_segment(
          ggplot2::aes(
            x = 0,
            xend = .data$logFC,
            yend = .data$gene_cell_key
          ),
          linewidth = 0.45,
          alpha = 0.70
        ) +
        ggplot2::geom_point(
          ggplot2::aes(
            size = .data$negative_log10_FDR
          ),
          alpha = 0.92
        ) +
        ggplot2::facet_wrap(
          ggplot2::vars(
            cell_type
          ),
          scales = "free_y",
          ncol = 3
        ) +
        ggplot2::scale_y_discrete(
          labels = function(x) {
            sub(
              "___.*$",
              "",
              x
            )
          }
        ) +
        ggplot2::scale_colour_manual(
          values = stats::setNames(
            c(
              "#2166AC",
              "#B2182B"
            ),
            c(
              paste0(
                "Higher in ",
                CONTROL_GROUP
              ),
              paste0(
                "Higher in ",
                CASE_GROUP
              )
            )
          ),
          name = NULL
        ) +
        ggplot2::scale_size_continuous(
          range = c(
            2,
            7
          ),
          name = "-log10(FDR)"
        ) +
        ggplot2::labs(
          title = "Statistically selected key genes",
          subtitle = paste0(
            "Top ",
            FINAL_EVIDENCE_TOP_GENES_PER_CELL_TYPE,
            " per cell type by FDR then absolute log2FC; display selection only"
          ),
          x = paste0(
            "log2 fold change (",
            CASE_GROUP,
            " versus ",
            CONTROL_GROUP,
            ")"
          ),
          y = NULL,
          caption = paste0(
            "Stress, mitochondrial and ribosomal genes are flagged in the ",
            "complete table but are not automatically removed."
          )
        ) +
        ggplot2::theme_minimal(
          base_size = 10
        ) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_blank(),
          panel.grid.major.y = ggplot2::element_line(
            colour = "#EEEEEE",
            linewidth = 0.25
          ),
          strip.text = ggplot2::element_text(
            face = "bold"
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          ),
          legend.position = "bottom"
        )

      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "final_key_genes_by_cell_type.png"
        ),
        p_final_key_genes,
        width = 14,
        height = max(
          9,
          3.2 *
            ceiling(
              dplyr::n_distinct(
                final_top_key_genes_plot$cell_type
              ) /
                3
            )
        ),
        dpi = 300
      )
    }


    # --------------------------------------------------------
    # 2.12 TARGETED PATHWAY-AXIS FIGURE
    # 2.12 定向通路轴图片
    # --------------------------------------------------------

    if (nrow(final_targeted_GSVA) > 0L) {

      final_targeted_plot_data <- final_targeted_GSVA |>
        dplyr::mutate(
          cell_type = factor(
            .data$cell_type,
            levels = final_cell_order
          ),
          axis_label = factor(
            .data$axis_label,
            levels = unname(
              FINAL_EVIDENCE_TARGETED_AXIS_LABELS[
                names(
                  FINAL_EVIDENCE_TARGETED_KEGG_AXES
                )
              ]
            )
          ),
          pathway_display = paste0(
            .data$pathway_name,
            " [",
            .data$pathway_id,
            "]"
          )
        )

      targeted_limit <- max(
        abs(
          final_targeted_plot_data$logFC
        ),
        na.rm = TRUE
      )

      if (
        !is.finite(
          targeted_limit
        ) ||
        targeted_limit == 0
      ) {
        targeted_limit <- 1
      }

      p_final_targeted_axes <- ggplot2::ggplot(
        final_targeted_plot_data,
        ggplot2::aes(
          x = .data$cell_type,
          y = .data$pathway_display,
          fill = .data$logFC
        )
      ) +
        ggplot2::geom_tile(
          colour = "white",
          linewidth = 0.35
        ) +
        ggplot2::geom_text(
          ggplot2::aes(
            label = ifelse(
              .data$significant,
              "*",
              ""
            )
          ),
          fontface = "bold",
          size = 4
        ) +
        ggplot2::facet_grid(
          rows = ggplot2::vars(
            axis_label
          ),
          scales = "free_y",
          space = "free_y"
        ) +
        ggplot2::scale_fill_gradient2(
          low = "#2166AC",
          mid = "white",
          high = "#B2182B",
          midpoint = 0,
          limits = c(
            -targeted_limit,
            targeted_limit
          ),
          na.value = "#BDBDBD",
          name = paste0(
            CASE_GROUP,
            " - ",
            CONTROL_GROUP,
            "\nGSVA logFC"
          )
        ) +
        ggplot2::labs(
          title = "Targeted KEGG pathway-axis activity",
          subtitle = paste0(
            "Secondary hypothesis-led view of predefined KEGG IDs; * = GSVA ",
            "FDR and effect thresholds both passed"
          ),
          x = NULL,
          y = NULL,
          caption = paste0(
            "A pathway association does not establish metabolite abundance, ",
            "disease diagnosis, malignant transformation or causal progression. ",
            "Grey = not tested or gene set unavailable."
          )
        ) +
        ggplot2::theme_minimal(
          base_size = 9
        ) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_text(
            angle = 45,
            hjust = 1
          ),
          strip.text.y = ggplot2::element_text(
            angle = 0,
            face = "bold",
            size = 8
          ),
          plot.title = ggplot2::element_text(
            face = "bold"
          )
        )

      save_plot_transparent(
        file.path(
          MDIR,
          "figures",
          "final_targeted_KEGG_axis_GSVA_heatmap.png"
        ),
        p_final_targeted_axes,
        width = max(
          13,
          0.65 *
            dplyr::n_distinct(
              final_targeted_plot_data$cell_type
            ) +
            5
        ),
        height = max(
          12,
          0.34 *
            dplyr::n_distinct(
              final_targeted_plot_data$pathway_display
            ) +
            6
        ),
        dpi = 300
      )
    }


    # --------------------------------------------------------
    # 2.13 DETERMINISTIC BILINGUAL REPORT
    # 2.13 规则驱动的双语报告
    # --------------------------------------------------------

    final_report_lines <- c(
      paste0(
        "# Final evidence synthesis / 最终证据整合：",
        DATASET_ID
      ),
      "",
      paste0(
        "Contrast / 比较：",
        CASE_GROUP,
        " versus ",
        CONTROL_GROUP
      ),
      "",
      paste0(
        "Generated / 生成时间：",
        format(
          Sys.time(),
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      "",
      "## Interpretation rules / 解释规则",
      "",
      paste0(
        "- Primary statistical families are pseudobulk DEG, pseudobulk GSVA ",
        "and the replicate-level exploratory composition screen. / 主要统计证据",
        "包括pseudobulk DEG、pseudobulk GSVA及重复层面的探索性细胞比例筛查。"
      ),
      paste0(
        "- GO and KEGG are interpretations of the same significant DEGs and are ",
        "not counted as independent confirmation. / GO和KEGG来源于同一批显著DEG，",
        "不作为独立重复验证。"
      ),
      paste0(
        "- Pseudotime and CellChat remain exploratory or descriptive. / 拟时与",
        "CellChat始终属于探索性或描述性结果。"
      ),
      paste0(
        "- SKIPPED or unavailable analyses are not evidence of no change. / ",
        "SKIPPED或不可用分析不代表不存在变化。"
      ),
      "",
      "## Cell-type evidence / Cell type证据"
    )

    for (
      i in seq_len(
        nrow(
          final_cell_evidence
        )
      )
    ) {

      row_i <- final_cell_evidence[i, ]

      top_gene_text <- final_key_genes |>
        dplyr::filter(
          .data$cell_type ==
            row_i$cell_type
        ) |>
        dplyr::arrange(
          .data$rank_within_cell_type
        ) |>
        dplyr::slice_head(
          n = 3L
        ) |>
        dplyr::pull(
          .data$gene
        ) |>
        paste(
          collapse = ", "
        )

      if (!nzchar(top_gene_text)) {
        top_gene_text <- "none / 无"
      }

      final_report_lines <- c(
        final_report_lines,
        "",
        paste0(
          "- **",
          row_i$cell_type,
          "** — ",
          row_i$evidence_level,
          "; significant DEG / 显著DEG = ",
          row_i$n_significant_DEGs,
          "; significant GSVA / 显著GSVA = ",
          row_i$n_significant_GSVA_pathways,
          "; composition FDR / 细胞比例FDR = ",
          ifelse(
            is.finite(
              row_i$composition_FDR
            ),
            formatC(
              row_i$composition_FDR,
              digits = 3,
              format = "g"
            ),
            "NA"
          ),
          "; pseudobulk = ",
          row_i$pseudobulk_status,
          "; top statistical genes / 前3个统计学基因 = ",
          top_gene_text,
          "."
        )
      )
    }

    final_skipped_cells <- final_cell_evidence |>
      dplyr::filter(
        .data$pseudobulk_status != "PASS"
      )

    final_report_lines <- c(
      final_report_lines,
      "",
      "## Explicit limitations / 明确局限性",
      "",
      paste0(
        "- Transcript abundance does not directly measure proteins, metabolites ",
        "or tissue histology. / 转录本丰度不能直接测量蛋白、代谢物或组织病理。"
      ),
      paste0(
        "- Targeted KEGG axes are hypothesis-led secondary summaries and do not ",
        "alter primary discovery. / 定向KEGG研究轴是预设假说驱动的二级汇总，",
        "不改变全局发现。"
      ),
      paste0(
        "- A cancer-related pathway name does not demonstrate malignant ",
        "transformation. / 癌症相关通路名称不能证明已经发生恶性转化。"
      )
    )

    if (nrow(final_skipped_cells) > 0L) {
      final_report_lines <- c(
        final_report_lines,
        "",
        "### Pseudobulk limitations / Pseudobulk限制"
      )

      for (
        i in seq_len(
          nrow(
            final_skipped_cells
          )
        )
      ) {
        final_report_lines <- c(
          final_report_lines,
          "",
          paste0(
            "- ",
            final_skipped_cells$cell_type[[i]],
            ": ",
            final_skipped_cells$pseudobulk_status[[i]],
            "; ",
            final_skipped_cells$pseudobulk_reason[[i]],
            ". No replicate-supported DEG conclusion is made. / 不作具有生物学",
            "重复支持的DEG结论。"
          )
        )
      }
    }

    final_report_lines <- c(
      final_report_lines,
      "",
      "## Output guide / 输出文件说明",
      "",
      paste0(
        "- `final_key_cell_types_evidence_matrix.csv`: complete evidence matrix ",
        "/ 完整cell type证据矩阵"
      ),
      paste0(
        "- `final_key_genes_statistical.csv`: all statistically selected genes ",
        "/ 全部统计学筛选基因"
      ),
      paste0(
        "- `final_key_genes_pathway_drivers.csv`: DEG-dependent pathway drivers ",
        "/ DEG依赖的通路驱动基因"
      ),
      paste0(
        "- `final_targeted_KEGG_axis_GSVA_results.csv`: hypothesis-led pathway ",
        "summary / 假说驱动的定向通路汇总"
      )
    )

    writeLines(
      final_report_lines,
      file.path(
        FINAL_REPORT_DIR,
        "final_evidence_summary_bilingual.md"
      )
    )
  }


  # ==========================================================
  # 3. FINAL OBJECT SUMMARY
  # ==========================================================

  final_summary <- data.frame(

    dataset_id =
      DATASET_ID,

    n_genes =
      nrow(obj),

    n_cells =
      ncol(obj),

    n_samples =
      length(
        unique(
          obj[[]][[SAMPLE_ID_COL]]
        )
      ),

    n_biological_replicates =
      length(
        unique(
          obj[[]][[BIOLOGICAL_REPLICATE_COL]]
        )
      ),

    n_conditions =
      length(
        unique(
          obj[[]][[CONDITION_COL]]
        )
      ),

    n_clusters =
      length(
        unique(
          obj$cluster
        )
      ),

    n_cell_types =
      length(
        unique(
          obj$cell_type
        )
      ),

    stringsAsFactors = FALSE
  )


  data.table::fwrite(
    final_summary,
    file.path(
      MDIR,
      "tables",
      "final_object_summary.csv"
    )
  )


  # ==========================================================
  # 4. SAVE FINAL RDS
  # ==========================================================

  final_file <- file.path(
    DIR_PROCESSED,
    paste0(
      DATASET_ID,
      "_FINAL_annotated_Seurat.rds"
    )
  )


  saveRDS(
    obj,
    final_file
  )


  # ==========================================================
  # 5. RESULT MANIFEST
  # ==========================================================

  result_manifest <- list.files(
    DIR_RESULTS,
    recursive = TRUE,
    full.names = TRUE
  )


  writeLines(
    result_manifest,
    file.path(
      MDIR,
      "logs",
      "result_manifest.txt"
    )
  )


  # ==========================================================
  # 6. SESSION INFORMATION
  # ==========================================================
  #
  # sessionInfo非常重要：
  # 可以记录R、Seurat、Bioconductor等package版本，
  # 对以后复现分析非常有价值。
  # ==========================================================

  writeLines(
    capture.output(
      sessionInfo()
    ),
    file.path(
      MDIR,
      "logs",
      "sessionInfo_end.txt"
    )
  )


  # ==========================================================
  # 7. FINAL MESSAGE
  # ==========================================================

  cat(
    "\n============================================================\n",
    "PIPELINE COMPLETED\n",
    "分析流程完成\n",
    "============================================================\n",
    "Dataset: ",
    DATASET_ID,
    "\n",
    "Genes: ",
    nrow(obj),
    "\n",
    "Cells: ",
    ncol(obj),
    "\n",
    "Clusters: ",
    length(unique(obj$cluster)),
    "\n",
    "Cell types: ",
    length(unique(obj$cell_type)),
    "\n\n",
    "Final object:\n",
    final_file,
    "\n",
    "============================================================\n",
    sep = ""
  )


  write_module_status(
    MODULE,
    "OK"
  )
}



# ============================================================
# RECOMMENDED EXECUTION AND REVIEW ORDER / 推荐运行与审阅顺序
# ============================================================
#
# Run the workflow in reviewable stages rather than executing M00-M15 blindly on
# a new dataset. Each stage has an explicit decision point:
# 新数据集不建议未经检查直接运行M00-M15，而应按可审阅阶段执行。每个阶段均设有
# 明确的继续条件：
#
# Stage 1 — M00-M03: environment, design and input provenance
# 阶段1——M00-M03：环境、实验设计与输入来源
#
#   Review package versions, sample-to-replicate mappings, group sizes, file
#   manifests and feature types. Continue only when sample identity and matrix
#   structure are unambiguous.
#   检查软件版本、样本—重复对应关系、各组重复数、文件清单及feature type。只有在
#   样本身份和矩阵结构均无歧义时才继续。
#
# Stage 2 — M04-M06: object construction and identity assignment
# 阶段2——M04-M06：对象构建与来源标签判定
#
#   Confirm RNA raw counts, optional HTO content, unique barcodes, sample metadata
#   and tissue assignments. Compare cell counts before and after any demultiplexing
#   or tissue subset.
#   确认RNA原始计数、可选HTO内容、barcode唯一性、样本metadata和组织标签。比较
#   demultiplexing或组织筛选前后的细胞数。
#
# Stage 3 — M07-M08: technical-quality filtering
# 阶段3——M07-M08：技术质量过滤
#
#   Inspect QC distributions and failure reasons by sample, then review predicted
#   doublet rates. Continue only when filtering is biologically plausible and no
#   replicate has been unintentionally depleted.
#   按样本检查QC分布和失败原因，再审阅预测doublet比例。只有在过滤符合生物学常识
#   且没有重复被意外大量删除时才继续。
#
# Stage 4 — M09-M10: exploratory representation and clustering
# 阶段4——M09-M10：探索性表达表征与聚类
#
#   Review normalization diagnostics, HVGs, PCA elbow/variance, cluster stability,
#   sample representation and marker coherence. Choose N_PCS_USE and resolution
#   from these diagnostics rather than from UMAP appearance alone.
#   检查标准化诊断、HVG、PCA elbow/方差、cluster稳定性、样本构成和marker一致性。
#   N_PCS_USE和resolution应依据这些诊断确定，而非仅凭UMAP外观。
#
# Stage 5 — M10B when justified: integration sensitivity analysis
# 阶段5——仅在有依据时运行M10B：整合敏感性分析
#
#   Enable only when there is evidence that technical sample effects obscure shared
#   identities. Compare integrated and non-integrated structures and verify that
#   disease-associated biology has not been erased.
#   仅在有证据表明技术性样本效应掩盖共同身份时启用。比较整合前后结构，并验证
#   疾病相关生物学信号未被消除。
#
# Stage 6 — M11: marker review and manual annotation
# 阶段6——M11：marker审阅与人工注释
#
#   Review top markers, canonical marker combinations, conflicting markers and
#   cluster quality. Complete manual_cluster_annotation.csv and rerun M11 so that
#   reviewed cell_type labels are written into the checkpoint.
#   检查top marker、经典marker组合、冲突marker和cluster质量。填写
#   manual_cluster_annotation.csv并重新运行M11，使审阅后的cell_type写入检查点。
#
# Stage 7 — M12-M13: replicate-level downstream analysis
# 阶段7——M12-M13：生物学重复层面的下游分析
#
#   First inspect descriptive composition by replicate. Then verify pseudobulk cell
#   and replicate minima, design-matrix rank and absence of perfect confounding
#   before interpreting edgeR results.
#   先按重复检查描述性细胞组成；随后在解释edgeR结果前，确认pseudobulk满足最少
#   细胞数和重复数、设计矩阵满秩且不存在完全混杂。
#
# Stage 8 — M14/M14B-M15: focused reporting, optional trajectory and handover
# 阶段8——M14/M14B-M15：重点报告、可选拟时与最终交付
#
#   Label post hoc focused analyses as exploratory. Run M14B only for a specified
#   plausible lineage, verify root/endpoint identity and representation across
#   biological replicates, and state that pseudotime is inferred rather than
#   observed time. Before handover, verify the final object, result manifest,
#   sample metadata, parameter record and software log together.
#   对事后选择的重点分析应明确标记为探索性。M14B只能用于明确指定且合理的谱系；
#   需要检查起点/终点身份及其在生物学重复中的覆盖，并说明拟时是推断顺序而非真实
#   时间。交付前需共同核对最终对象、结果清单、样本metadata、参数记录和软件日志。
#
# ============================================================


# ============================================================
# QUESTIONS AFTER EVERY MODULE / 每个模块结束后的审阅问题
# ============================================================
#
# 1. What exact object or table entered the module, and from which checkpoint?
#    本模块的精确输入对象或表格是什么，来自哪个检查点？
#
# 2. Which data representation was used: raw counts, normalized expression,
#    scaled values or an integrated representation? Was that representation
#    appropriate for the stated purpose?
#    本模块使用的是原始计数、标准化表达、scaled值还是整合表示？该数据表示是否适合
#    当前目的？
#
# 3. What changed in the object—cells, genes, assays/layers, metadata, identities
#    or reductions—and was the change recorded in an output table or log?
#    对象中的细胞、基因、assay/layer、metadata、identity或降维结果发生了什么变化？
#    这些变化是否记录在输出表或日志中？
#
# 4. Which parameters affected the result, which values are generic defaults and
#    which were justified from this dataset's diagnostics?
#    哪些参数会影响结果？哪些是通用默认值，哪些已经由本数据集的诊断结果支持？
#
# 5. Is the output a technical diagnostic, an exploratory biological description
#    or a formal statistical inference? Is it labelled accordingly?
#    当前输出属于技术诊断、探索性生物学描述，还是正式统计推断？报告中的表述是否与
#    其证据等级一致？
#
# 6. What is the independent unit at this step: cell, technical library or
#    biological replicate? Could pseudoreplication affect the conclusion?
#    本步骤的独立单位是细胞、技术文库还是生物学重复？结论是否可能受到伪重复影响？
#
# 7. Are apparent condition differences potentially explained by sample imbalance,
#    batch, cell recovery, QC, doublet removal or annotation uncertainty?
#    表面condition差异是否可能由样本不平衡、batch、细胞捕获量、QC、doublet删除或
#    注释不确定性解释？
#
# 8. Did every expected replicate remain represented, and is any result dominated
#    by a single replicate or unusually large library?
#    所有预期生物学重复是否仍然存在？结果是否被单个重复或异常大的文库主导？
#
# 9. If integration was used, did it preserve condition-associated states and
#    remain excluded from count-based statistical testing?
#    若使用整合，是否保留了condition相关状态，并且未将整合值用于基于计数的统计检验？
#
# 10. Are the module's principal tables, figures, warnings and status log present,
#     and can the next module be rerun from the saved checkpoint without relying on
#     undocumented objects in the current R session?
#     本模块的主要表格、图片、warning和状态日志是否齐全？下一模块能否仅依赖已保存的
#     检查点重新运行，而不依赖当前R会话中未记录的对象？
#
# ============================================================


# ============================================================
# END
# ============================================================
