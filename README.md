# SELEX DNA 文库序列分析平台

## 项目简介

**SELEX DNA 文库序列分析平台** 是一个面向 SELEX DNA paired-end 测序数据的 Web 分析平台。

平台将原有的序列处理、文库归类、序列计数、CD-HIT-EST 家族分析和 DNA motif 可视化流程整合到一个网页中，使用户可以通过浏览器完成数据上传、参数设置、分析运行和结果下载。

---

## 一、主要功能

### 1. 测序数据处理

支持上传 paired-end 测序数据：

- R1 FASTQ / FASTQ.GZ
- R2 FASTQ / FASTQ.GZ

平台按照既定分析流程进行：

**测序数据 → 去接头 → 双端序列合并 → 引物识别 → 随机区提取**

### 2. Tag / 无 Tag 文库处理

平台支持有 Tag 和无 Tag 两种文库情况。

用户可以在网页中填写 Tag 信息：

- 有 Tag：根据输入的 Tag 信息对序列进行文库归类
- 无 Tag：Tag 信息留空，序列进入无 Tag 文库分析

其中：

> `random_length` 为用户手动设定的随机序列长度，与是否存在 Tag 无关。

### 3. 序列计数与过滤

对提取后的 DNA 序列进行计数，并按照设定的 counts 阈值进行过滤，减少低丰度序列对后续分析的影响。

### 4. DNA 家族分析

使用 **CD-HIT-EST** 对 DNA 序列进行家族聚类。

当前网页默认示例参数：

```text
CD-HIT-EST 相似度阈值：80%
```

通过家族聚类可以将具有较高序列相似性的 DNA 序列归入同一家族，并进一步统计各家族的成员数量。

### 5. DNA Motif 分析

根据家族成员数量对 DNA 家族进行排序，并对高丰度家族进行多序列比对和 motif 分析。

最终生成：

- DNA motif letter plot
- PNG 图像
- PDF 图像
- 家族统计结果
- 代表序列及相关分析结果

---

# 二、网站分析流程

```text
        ┌─────────────────────┐
        │   上传 R1 / R2 数据  │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │     去除测序接头     │
        │       Cutadapt      │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │   paired-end 合并    │
        │        FLASH        │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │   引物识别 / 提取    │
        │   前引物 + 后引物     │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │    提取随机序列区    │
        │   random_length      │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │   Tag / 无 Tag 归类   │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │     counts 过滤     │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │   CD-HIT-EST 聚类    │
        │      80% 相似度      │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │   Family 数量排序    │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │   DNA Motif 分析     │
        │    Letter Plot      │
        └──────────┬──────────┘
                   ↓
        ┌─────────────────────┐
        │      结果下载        │
        └─────────────────────┘
```

---

# 三、网站主要输入参数

| 参数 | 作用 |
|---|---|
| R1 文件 | paired-end 测序正向数据 |
| R2 文件 | paired-end 测序反向数据 |
| Adapter 1 | R1 接头序列 |
| Adapter 2 | R2 接头序列 |
| Forward primer | 前引物 |
| Reverse primer | 后引物 |
| Random length | 随机序列长度 |
| Tag | 文库 Tag 信息 |
| PCR reverse primer | 可选 PCR reverse primer |
| FLASH overlap | 双端序列拼接参数 |
| CD-HIT threshold | DNA 家族相似度阈值 |
| Threads | 分析线程数 |
| Min counts | 序列最低 counts |
| Family member threshold | 进行 motif 分析的家族成员数阈值 |
| Single-sequence count threshold | 单序列丰度阈值 |
| Motif coverage | motif 分析覆盖度参数 |

---

# 四、项目结构

```text
SELEX_GitHub_Complete/
│
├── app.py
│   └── FastAPI 网站后端/API
│
├── analysis_engine.py
│   └── SELEX DNA 分析流程
│
├── analysis_engine_original.py
│   └── 原始分析代码备份
│
├── static/
│   ├── index.html
│   │   └── 网站首页
│   ├── style.css
│   │   └── 网站样式
│   └── script.js
│       └── 网页交互、任务状态和结果下载
│
├── Dockerfile
│   └── 网站运行环境
│
├── requirements.txt
│   └── Python 依赖
│
├── render.yaml
│   └── Render 部署配置
│
├── jobs/
│   └── 分析任务临时目录
│
├── data/
│   └── 数据目录
│
└── README.md
    └── 项目说明
```

---

# 五、技术组成

本平台主要由以下部分组成：

```text
前端
HTML + CSS + JavaScript
        ↓
后端
FastAPI
        ↓
分析程序
Python
        ↓
生物信息学工具
Cutadapt
FLASH
CD-HIT-EST
        ↓
结果
Family analysis
DNA Motif
Letter Plot
```

同时使用：

- Biopython
- regex
- pandas
- NumPy
- Matplotlib
- Logomaker

---

# 六、运行方式

## 方式 1：在线部署

项目可以部署到支持 Docker 的云平台，例如 Render。

部署完成后，用户通过浏览器访问网站即可使用分析平台。

网站本身负责：

1. 上传测序数据
2. 提交分析参数
3. 启动分析任务
4. 查看分析进度
5. 下载分析结果

## 方式 2：本地 Docker 运行

在项目根目录执行：

```bash
docker build -t selex-web .
docker run --rm -p 10000:10000 selex-web
```

然后浏览器访问：

```text
http://localhost:10000
```

---

# 七、数据安全说明

真实测序数据不应直接提交到 GitHub。

项目的 `.gitignore` 已经忽略常见测序文件：

```text
*.fastq
*.fastq.gz
*.fq
*.fq.gz
```

真实 FASTQ 数据应通过网站上传或在本地运行分析程序。

---

# 八、当前版本说明

当前版本主要用于：

- 课题组内部展示
- SELEX 分析流程验证
- 小规模数据测试
- Web 平台原型部署

当前代码中，Tag / 无 Tag 模式根据网页提交的 Tag 列表进行处理；`random_length` 始终由用户手动设定，不根据 Tag 状态自动改变。

对于大规模数据、长时间运行任务或多人同时使用的正式生产环境，后续可进一步增加：

- 后台任务队列
- 独立分析 Worker
- 持久化结果存储
- 分析历史记录
- 用户登录与权限管理
- 结果 ZIP 一键下载
- Family counts 可视化
- Top family DNA motif 网页预览

---

## 项目目标

最终形成一个适用于课题组 SELEX DNA 数据分析的统一 Web 平台：

**数据上传 → 参数设置 → 自动分析 → Family 分析 → DNA Motif 可视化 → 结果下载**
