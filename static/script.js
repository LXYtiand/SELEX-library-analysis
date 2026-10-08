const form = document.getElementById("analysisForm");
const result = document.getElementById("result");
const statusBanner = document.getElementById("statusBanner");
const jobSummary = document.getElementById("jobSummary");
const logBox = document.getElementById("log");
const filesBox = document.getElementById("files");
const progressBar = document.getElementById("progressBar");
const startButton = document.getElementById("startButton");
const tagsBox = document.getElementById("tags");
const tagStatus = document.getElementById("tagStatus");
const tagCount = document.getElementById("tagCount");

const stages = [
  ["cutadapt", ["Cutadapt", "去接头"]],
  ["flash", ["FLASH", "拼接"]],
  ["extract", ["随机区", "序列提取", "文库"]],
  ["family", ["CD-HIT", "家族", "family"]],
  ["motif", ["motif", "logo", "字母图"]]
];

document.querySelectorAll('input[type="file"]').forEach(input => {
  input.addEventListener("change", () => {
    const target = document.querySelector(`[data-for="${input.name}"]`);
    if (target) target.textContent = input.files?.[0]?.name || "未选择文件";
  });
});

function parseTags() {
  return tagsBox.value
    .split(/[\n,，]+/)
    .map(x => x.trim())
    .filter(Boolean);
}

function updateTagStatus() {
  const count = parseTags().length;
  tagCount.textContent = `${count} 个 Tag`;
  if (count > 0) {
    tagStatus.innerHTML = '<span class="status-dot"></span>当前：有 Tag · 将按文库字段归类';
    tagStatus.style.color = "var(--success)";
  } else {
    tagStatus.innerHTML = '<span class="status-dot"></span>当前：无 Tag · 所有序列进入 library_1';
    tagStatus.style.color = "var(--muted)";
  }
}
tagsBox.addEventListener("input", updateTagStatus);
updateTagStatus();

function setStage(name, state) {
  const el = document.querySelector(`.stage[data-stage="${name}"]`);
  if (!el) return;
  el.classList.remove("active", "done");
  if (state) el.classList.add(state);
}

function updateStages(log) {
  const text = (log || "").toLowerCase();
  const checks = [
    ["cutadapt", ["cutadapt", "去接头", "过滤与去接头"]],
    ["flash", ["flash", "拼接双端"]],
    ["extract", ["提取", "随机区", "library_"]],
    ["family", ["cd-hit", "cdhit", "家族", "family_statistics"]],
    ["motif", ["motif", "logo", "sequence logo"]]
  ];
  let latest = -1;
  checks.forEach((item, idx) => {
    const found = item[1].some(k => text.includes(k.toLowerCase()));
    if (found) latest = Math.max(latest, idx);
  });

  checks.forEach((item, idx) => {
    if (idx < latest) setStage(item[0], "done");
    else if (idx === latest) setStage(item[0], "active");
    else setStage(item[0], "");
  });

  const progress = latest < 0 ? 3 : Math.min(95, 10 + latest * 21);
  progressBar.style.width = `${progress}%`;
}

async function getJSON(url, options) {
  const res = await fetch(url, options);
  let data = {};
  try { data = await res.json(); } catch (_) {}
  if (!res.ok) throw new Error(data.detail || data.error || `请求失败 (${res.status})`);
  return data;
}

function showStatus(text, type="") {
  statusBanner.textContent = text;
  statusBanner.className = "status-banner" + (type ? ` ${type}` : "");
}

form.addEventListener("submit", async (e) => {
  e.preventDefault();

  result.classList.remove("hidden");
  filesBox.innerHTML = "";
  logBox.textContent = "";
  progressBar.style.width = "3%";
  document.querySelectorAll(".stage").forEach(x => x.classList.remove("active","done"));
  showStatus("正在上传测序数据并创建分析任务……");
  jobSummary.textContent = "任务启动中，请稍候。";
  startButton.disabled = true;

  try {
    const fd = new FormData(form);
    const data = await getJSON("/api/analyze", { method: "POST", body: fd });

    const modeText = data.mode || "未说明";
    jobSummary.textContent = `任务 ID：${data.job_id} · 模式：${modeText}`;
    showStatus(`任务已启动：${data.job_id} · ${modeText}`);
    await poll(data.job_id);
  } catch (err) {
    showStatus(`启动失败：${err.message}`, "error");
    jobSummary.textContent = "任务未启动";
    startButton.disabled = false;
  }
});

async function poll(jobId) {
  try {
    const status = await getJSON(`/api/jobs/${jobId}`);
    const logData = await fetch(`/api/jobs/${jobId}/log`)
      .then(r => r.ok ? r.json() : {log:""})
      .catch(() => ({log:""}));

    logBox.textContent = logData.log || "";
    logBox.scrollTop = logBox.scrollHeight;
    updateStages(logData.log || "");

    if (status.status === "running" || status.status === "queued") {
      showStatus(status.status === "queued" ? "任务排队中……" : "分析正在运行……");
      setTimeout(() => poll(jobId), 3000);
      return;
    }

    if (status.status === "completed") {
      progressBar.style.width = "100%";
      document.querySelectorAll(".stage").forEach(x => x.classList.add("done"));
      showStatus("分析完成，可以下载结果。", "success");
      jobSummary.textContent = `任务 ID：${jobId} · 分析完成`;
    } else {
      progressBar.style.width = "100%";
      showStatus(`分析失败：${status.error || "请查看运行日志"}`, "error");
      jobSummary.textContent = `任务 ID：${jobId} · 分析失败`;
    }

    const files = await getJSON(`/api/jobs/${jobId}/files`);
    renderFiles(jobId, files.files || []);
  } catch (err) {
    showStatus(`状态查询失败：${err.message}`, "error");
  } finally {
    startButton.disabled = false;
  }
}

function renderFiles(jobId, files) {
  if (!files.length) {
    filesBox.innerHTML = "<h3>结果文件</h3><p class='file-name'>暂未发现可下载的结果文件，请查看运行日志。</p>";
    return;
  }

  const rows = files.map(path => {
    const href = `/api/jobs/${encodeURIComponent(jobId)}/download/${path.split("/").map(encodeURIComponent).join("/")}`;
    const name = path.split("/").pop();
    return `<div class="file"><span>📄</span><a href="${href}" target="_blank" rel="noopener">${path}</a></div>`;
  }).join("");

  filesBox.innerHTML = `<h3>结果文件（${files.length}）</h3>${rows}`;
}
