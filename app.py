import json
import os
import shutil
import subprocess
import sys
import uuid
from pathlib import Path
from threading import Thread

from fastapi import FastAPI, File, Form, UploadFile, HTTPException
from fastapi.responses import FileResponse, HTMLResponse
from fastapi.staticfiles import StaticFiles

BASE_DIR = Path(__file__).resolve().parent
JOBS_DIR = Path(os.environ.get("JOBS_DIR", BASE_DIR / "jobs"))
JOBS_DIR.mkdir(parents=True, exist_ok=True)

app = FastAPI(title="SELEX Library Sequence Analysis")

app.mount("/static", StaticFiles(directory=BASE_DIR / "static"), name="static")


def save_upload(upload: UploadFile, path: Path):
    with path.open("wb") as f:
        while True:
            chunk = upload.file.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)


def run_job(job_dir: Path):
    config_path = job_dir / "config.json"
    log_path = job_dir / "run.log"
    status_path = job_dir / "status.json"

    status_path.write_text(json.dumps({"status": "running"}, ensure_ascii=False), encoding="utf-8")

    env = os.environ.copy()
    env["SELEX_CONFIG_PATH"] = str(config_path)

    with log_path.open("w", encoding="utf-8") as log:
        try:
            result = subprocess.run(
                [sys.executable, str(BASE_DIR / "analysis_engine.py")],
                cwd=str(job_dir),
                env=env,
                stdout=log,
                stderr=subprocess.STDOUT,
                text=True,
            )
            status = "completed" if result.returncode == 0 else "failed"
            status_path.write_text(
                json.dumps(
                    {"status": status, "returncode": result.returncode},
                    ensure_ascii=False
                ),
                encoding="utf-8",
            )
        except Exception as exc:
            with log.open("a", encoding="utf-8") as log2:
                log2.write(f"\nSERVER ERROR: {exc}\n")
            status_path.write_text(
                json.dumps({"status": "failed", "error": str(exc)}, ensure_ascii=False),
                encoding="utf-8",
            )


@app.get("/", response_class=HTMLResponse)
def index():
    return FileResponse(BASE_DIR / "static" / "index.html")


@app.get("/health")
def health():
    return {"status": "ok"}


@app.post("/api/analyze")
async def analyze(
    r1: UploadFile = File(...),
    r2: UploadFile = File(...),
    adapter1: str = Form(...),
    adapter2: str = Form(...),
    extract_begin_seq: str = Form(...),
    extract_end_seq: str = Form(...),
    random_length: int = Form(...),
    tags: str = Form(""),
    reverse_primer: str = Form(""),
    remove_pcr_primer: str = Form("on"),
    max_mismatch: int = Form(2),
    flash_max_overlap: int = Form(250),
    sim_threshold: float = Form(0.80),
    threads: int = Form(8),
    min_retention_count: int = Form(3),
    member_count_thresh: int = Form(50),
    single_seq_count_thresh: int = Form(500),
    min_col_coverage: float = Form(0.30),
):
    if random_length <= 0:
        raise HTTPException(400, "random_length 必须大于 0")
    if not (0 < sim_threshold <= 1):
        raise HTTPException(400, "CD-HIT 相似度阈值必须在 0~1 之间")
    if max_mismatch < 0:
        raise HTTPException(400, "最大错配数不能小于 0")
    if threads < 1:
        raise HTTPException(400, "线程数必须大于等于 1")
    if min_retention_count < 1:
        raise HTTPException(400, "最小 counts 必须大于等于 1")
    if not (0 <= min_col_coverage <= 1):
        raise HTTPException(400, "最小位点覆盖度必须在 0~1 之间")

    job_id = uuid.uuid4().hex[:12]
    job_dir = JOBS_DIR / job_id
    job_dir.mkdir(parents=True)

    r1_path = job_dir / "R1.fastq.gz"
    r2_path = job_dir / "R2.fastq.gz"
    save_upload(r1, r1_path)
    save_upload(r2, r2_path)

    tag_list = [x.strip().upper() for x in tags.replace(",", "\n").splitlines() if x.strip()]

    config = {
        "input_file1": str(r1_path),
        "input_file2": str(r2_path),
        "final_output_dir": str(job_dir / "output_detail"),
        "OUTPUT_DIR": str(job_dir / "family_results"),
        "output_prefix": str(job_dir / "processed"),
        "adapter1": adapter1.strip().upper(),
        "adapter2": adapter2.strip().upper(),
        "extract_begin_seq": extract_begin_seq.strip().upper(),
        "extract_end_seq": extract_end_seq.strip().upper(),
        "random_length": random_length,
        "second_half_fields": tag_list,
        "reverse_primer": reverse_primer.strip().upper(),
        "REMOVE_PCR_PRIMER": str(remove_pcr_primer).lower() in {"on", "true", "1", "yes"},
        "max_mismatch": max_mismatch,
        "FLASH_MAX_OVERLAP": str(flash_max_overlap),
        "SIM_THRESHOLD": sim_threshold,
        "THREADS": threads,
        "MIN_RETENTION_COUNT": min_retention_count,
        "MEMBER_COUNT_THRESH": member_count_thresh,
        "SINGLE_SEQ_COUNT_THRESH": single_seq_count_thresh,
        "MIN_COL_COVERAGE": min_col_coverage,
        "select_file_names": None,
    }
    (job_dir / "config.json").write_text(
        json.dumps(config, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    Thread(target=run_job, args=(job_dir,), daemon=True).start()

    return {
        "job_id": job_id,
        "mode": "有 Tag" if tag_list else "无 Tag",
        "tag_count": len(tag_list),
        "status": "running",
    }


@app.get("/api/jobs/{job_id}")
def job_status(job_id: str):
    job_dir = JOBS_DIR / job_id
    status_path = job_dir / "status.json"
    if not job_dir.exists():
        raise HTTPException(404, "找不到任务")
    if status_path.exists():
        status = json.loads(status_path.read_text(encoding="utf-8"))
    else:
        status = {"status": "queued"}
    return status


@app.get("/api/jobs/{job_id}/log")
def job_log(job_id: str):
    path = JOBS_DIR / job_id / "run.log"
    if not path.exists():
        raise HTTPException(404, "日志尚未生成")
    return {"log": path.read_text(encoding="utf-8", errors="ignore")[-30000:]}


@app.get("/api/jobs/{job_id}/files")
def job_files(job_id: str):
    job_dir = JOBS_DIR / job_id
    if not job_dir.exists():
        raise HTTPException(404, "找不到任务")
    files = []
    for p in job_dir.rglob("*"):
        if p.is_file() and p.name not in {"config.json", "status.json", "run.log", "R1.fastq.gz", "R2.fastq.gz"}:
            files.append(str(p.relative_to(job_dir)))
    return {"files": files}


@app.get("/api/jobs/{job_id}/download/{file_path:path}")
def download_file(job_id: str, file_path: str):
    job_dir = (JOBS_DIR / job_id).resolve()
    target = (job_dir / file_path).resolve()
    if job_dir not in target.parents or not target.is_file():
        raise HTTPException(404, "文件不存在")
    return FileResponse(target)


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host="0.0.0.0", port=int(os.environ.get("PORT", "10000")))
