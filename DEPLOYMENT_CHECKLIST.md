# 部署检查清单

## 本地

- [ ] ZIP 已解压
- [ ] `index.html`、`style.css`、`script.js` 存在
- [ ] Docker 已安装
- [ ] `docker build -t selex-web .` 成功
- [ ] `docker run --rm -p 10000:10000 selex-web` 成功
- [ ] 浏览器可以打开 `http://localhost:10000`
- [ ] `/health` 返回 `{"status":"ok"}`

## GitHub

- [ ] 已创建 `SELEX-library-analysis`
- [ ] 没有提交 FASTQ/FASTQ.GZ
- [ ] 已执行 `git add .`
- [ ] 已执行 `git commit`
- [ ] 已执行 `git push -u origin main`

## Render

- [ ] GitHub 已连接 Render
- [ ] Runtime 使用 Docker
- [ ] Dockerfile 为项目根目录的 `Dockerfile`
- [ ] Health Check Path 为 `/health`
- [ ] 部署成功
- [ ] 在线网址可以打开
- [ ] 用小型测试 FASTQ 做一次完整流程测试
