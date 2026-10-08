# SELEX DNA 文库序列分析平台

这是一个针对 **SELEX DNA paired-end 测序数据** 的 Web 部署项目。

## 1. 项目结构

```text
SELEX_GitHub_Complete/
├── app.py                         # FastAPI 网页后端/API
├── analysis_engine.py             # 原 SELEX 分析流程（网页配置版）
├── analysis_engine_original.py    # 原始分析代码备份
├── static/
│   ├── index.html                 # 网站页面
│   ├── style.css                  # 网站样式
│   └── script.js                  # 上传、任务状态、结果下载
├── Dockerfile                     # Docker 运行环境
├── requirements.txt               # Python 依赖
├── render.yaml                    # Render 配置
├── README.md
├── .gitignore
├── jobs/.gitkeep                  # 临时任务目录
└── data/.gitkeep
```

## 2. 当前实际分析流程

网站对应现有 Python 分析程序：

1. 上传 R1/R2 FASTQ 或 FASTQ.GZ
2. Cutadapt 去接头
3. FLASH 拼接 paired-end reads
4. 可选去除 PCR reverse primer
5. 根据用户手动设定的 `random_length` 提取随机区
6. Tag 填写时：按 Tag 对应文库归类；Tag 留空时：无 Tag，统一进入 `library_1`
7. counts 过滤
8. CD-HIT-EST 家族聚类，默认相似度 80%
9. 输出家族统计、代表序列、按丰度排序的序列
10. 多序列比对并生成 DNA motif logo PNG/PDF

**注意：** 当前代码的 Tag/无 Tag 判断依据是网页是否提交 Tag 列表；`random_length` 始终独立于 Tag 模式。

## 3. 本地测试

推荐先安装 Docker。

在项目根目录执行：

```bash
docker build -t selex-web .
docker run --rm -p 10000:10000 selex-web
```

浏览器打开：

```text
http://localhost:10000
```

如果需要不使用 Docker 的本地测试，也可以安装 requirements 后运行：

```bash
pip install -r requirements.txt
uvicorn app:app --host 127.0.0.1 --port 10000
```

注意：本地非 Docker 运行还必须自行安装并确保以下命令在 PATH 中：

```text
cutadapt
flash
cd-hit-est
```

## 4. 上传到 GitHub

### 第一步：网页创建 Repository

GitHub → New repository

推荐名称：

```text
SELEX-library-analysis
```

创建时不要再自动生成 README/.gitignore，因为项目里已经有。

### 第二步：PowerShell 进入项目目录

```powershell
cd D:\你的路径\SELEX_GitHub_Complete
```

### 第三步：初始化并上传

```powershell
git init
git add .
git commit -m "Initial SELEX web platform"
git branch -M main
git remote add origin https://github.com/YOUR_USERNAME/SELEX-library-analysis.git
git push -u origin main
```

## 5. Render 部署

Render 不需要下载安装，直接使用网页版。

进入 Render → New → Web Service → 连接 GitHub Repository。

本项目使用 Dockerfile：

```text
Dockerfile
```

Health Check：

```text
/health
```

部署完成后 Render 会提供一个在线网址。

## 6. 非常重要：不要上传真实 FASTQ

`.gitignore` 已经忽略：

```text
*.fastq
*.fastq.gz
*.fq
*.fq.gz
```

真实测序数据请通过网站上传，不要提交到 GitHub。

## 7. 当前版本的部署边界

这是一个可以用于演示和小规模测试的完整 Web 原型。长时间、大文件、高并发的正式课题组使用，建议后续升级为：

```text
Web/API
   ↓
任务队列
   ↓
独立分析 Worker
   ↓
持久化磁盘/对象存储
```

因为 Render 普通 Web Service 的本地任务目录不应作为长期结果存储。

## 8. 推荐的后续开发

可以继续增加：

- 自动展示总 reads、有效序列数、library 数
- Family counts 排名表
- Top family DNA motif logo 网页预览
- 一键下载全部结果 ZIP
- 分析历史记录
- 用户登录与权限
- 持久化结果存储
