#!/usr/bin/env python3
"""
分子对接引擎 — Vina 对接 + 3D 可视化数据
参考：AutoDock Vina 官方文档，默认 exhaustiveness=8
"""
import sqlite3, os, subprocess, urllib.request, traceback, datetime

BASE = os.path.dirname(os.path.abspath(__file__))
LOG_PATH = os.path.join(BASE, "_dock_log.txt")

def _log(msg):
    """写日志到文件，绕过 Flask stdout 缓冲问题"""
    try:
        with open(LOG_PATH, "a", encoding="utf-8") as f:
            f.write(f"{datetime.datetime.now().strftime('%H:%M:%S')} {msg}\n")
    except: pass

# ── 工具路径 ─────────────────────────────────────
OBABEL = r"C:\Program Files\OpenBabel-2.4.1\obabel.exe"
if not os.path.exists(OBABEL):
    for alt in [r"C:\Program Files\OpenBabel\obabel.exe", r"C:\OpenBabel\obabel.exe"]:
        if os.path.exists(alt): OBABEL = alt; break

VINA = None
for p in [os.path.join(BASE, "vina.exe"), r"D:\autodock\vina.exe"]:
    if os.path.exists(p): VINA = p; break
if not VINA and os.path.exists(r"D:\autodock"):
    files = sorted([f for f in os.listdir(r"D:\autodock")
                    if os.path.isfile(os.path.join(r"D:\autodock", f))])
    for f in files:
        fl = f.lower()
        if fl in ('vina.exe', 'vina.exe.exe'):
            VINA = os.path.join(r"D:\autodock", f); break
    if not VINA:
        for f in files:
            fl = f.lower()
            if fl.endswith('.exe') and 'vina' in fl and 'split' not in fl:
                VINA = os.path.join(r"D:\autodock", f); break

# 预置口袋中心（基于 PDB 共晶结构人工标定）
# 口袋中心（仅保留手动验证过的）
POCKET_CTR = {
    "3O96": (22.5, 8.5, 22.0),  # AKT1 ATP口袋
}


def _download_pdb(pdb_id, save_path):
    import time as _time
    for attempt in range(3):
        try:
            req = urllib.request.Request(
                f"https://files.rcsb.org/download/{pdb_id}.pdb",
                headers={'User-Agent': 'Mozilla/5.0'}
            )
            with urllib.request.urlopen(req, timeout=60) as r:
                data = r.read()
                if len(data) > 3000:
                    with open(save_path, "wb") as f: f.write(data)
                    return True
        except Exception as e:
            if attempt < 2:
                _log(f"RCSB {pdb_id} attempt {attempt+1} failed, retrying...")
                _time.sleep(2)
            else:
                _log(f"RCSB {pdb_id} failed after 3 attempts: {e}")
    return False


def _download_alphafold(uniprot, save_path):
    for ver in ["v4", "v3"]:
        try:
            req = urllib.request.Request(
                f"https://alphafold.ebi.ac.uk/files/AF-{uniprot}-F1-model_{ver}.pdb",
                headers={'User-Agent': 'Mozilla/5.0'}
            )
            with urllib.request.urlopen(req, timeout=120) as r:
                data = r.read()
                if len(data) > 3000:
                    with open(save_path, "wb") as f: f.write(data)
                    _log(f"AlphaFold {ver} OK for {uniprot}")
                    return True
        except Exception as e:
            _log(f"AlphaFold {ver} failed: {e}")
    return False


def _receptor_fragment(pdb_path, center, radius=30):
    """提取口袋附近原子用于 3D 可视化"""
    cx, cy, cz = center; r2 = radius * radius; kept = []
    try:
        with open(pdb_path) as f:
            for line in f:
                if line.startswith("ATOM") or line.startswith("HETATM"):
                    try:
                        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        if (x-cx)**2 + (y-cy)**2 + (z-cz)**2 < r2: kept.append(line)
                    except: pass
                elif line.startswith("TER") or line.startswith("END"): kept.append(line)
    except: pass
    return "".join(kept[:50000])


def dock_fast(smiles, uniprot, receptor_pdb=None, force_blind=False, box_size=28):
    _log(f"START {uniprot} {smiles[:50]}")

    result = {
        "status": "error", "gene": "", "affinity": None, "method": "unknown",
        "pdb": "", "pose_sdf": "", "receptor_pdb": "", "pocket_center": [0, 0, 0],
        "ligand_smiles": smiles,
    }

    db_path = os.path.join(BASE, "tcm_demo.db")
    if not os.path.exists(db_path): return result

    conn = sqlite3.connect(db_path); conn.row_factory = sqlite3.Row; c = conn.cursor()
    c.execute("SELECT gene_name, pocket_score, has_pocket FROM targets WHERE uniprot_id=?", (uniprot,))
    row = c.fetchone(); conn.close()
    if not row: return result

    gene = row["gene_name"]; ps = row["pocket_score"] if row["has_pocket"] else 0
    result["gene"] = gene
    fpocket = round(-6.0 - ps * 4.0, 2) if ps >= 0.3 else (round(-4.0 - ps * 3.0, 2) if ps > 0 else None)

    if not VINA or not os.path.exists(OBABEL):
        result.update({"status": "ok", "affinity": fpocket, "method": "Fpocket (no tools)"})
        return result

    # ── 1. 配体 PDBQT ────────────────────────────
    import hashlib
    _uid = hashlib.md5((smiles+uniprot+datetime.datetime.now().isoformat()).encode()).hexdigest()[:8]
    lig = os.path.join(BASE, f"_dl_{_uid}.pdbqt")
    smi_file = os.path.join(BASE, f"_dl_{_uid}.smi")
    mol2_file = os.path.join(BASE, f"_dl_{_uid}.mol2")
    with open(smi_file, "w") as f: f.write(smiles)
    success = False
    # 1: 标准直接转（最稳定）
    try:
        subprocess.run([OBABEL, smi_file, "-opdbqt", "--gen3d", "-O", lig], timeout=180, capture_output=True)
        if os.path.exists(lig) and os.path.getsize(lig) >= 100: success = True
    except: pass
    # 2: 不用gen3d，obabel自动生成
    if not success:
        try:
            subprocess.run([OBABEL, smi_file, "-opdbqt", "-O", lig], timeout=120, capture_output=True)
            if os.path.exists(lig) and os.path.getsize(lig) >= 100: success = True
        except: pass
    # 3: 两步法 SMILES→MOL2→PDBQT
    if not success:
        try:
            subprocess.run([OBABEL, smi_file, "-omol2", "--gen3d", "-O", mol2_file], timeout=180, capture_output=True)
            if os.path.exists(mol2_file) and os.path.getsize(mol2_file) >= 100:
                subprocess.run([OBABEL, mol2_file, "-opdbqt", "-O", lig], timeout=60, capture_output=True)
                if os.path.exists(lig) and os.path.getsize(lig) >= 100: success = True
        except: pass
    try: os.remove(smi_file)
    except: pass
    try: os.remove(mol2_file)
    except: pass
    # 4: RDKit ETKDGv3（复杂天然产物）
    if not success:
        _log(f"Ligand: trying RDKit for {uniprot}")
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                mol = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 42
                cid = AllChem.EmbedMolecule(mol, params)
                _log(f"RDKit embed result={cid} natoms={mol.GetNumAtoms()}")
                if cid >= 0:
                    mol = Chem.RemoveHs(mol)
                    rd_sdf = os.path.join(BASE, f"_dl_{_uid}_r.sdf")
                    writer = Chem.SDWriter(rd_sdf)
                    writer.write(mol)
                    writer.close()
                    subprocess.run([OBABEL, rd_sdf, "-opdbqt", "-O", lig], timeout=60, capture_output=True)
                    if os.path.exists(lig) and os.path.getsize(lig) >= 100:
                        success = True
                        _log(f"RDKit ligand OK: {os.path.getsize(lig)}B")
                    else:
                        _log(f"RDKit→PDBQT failed: size={os.path.getsize(lig) if os.path.exists(lig) else 0}")
                    try: os.remove(rd_sdf)
                    except: pass
        except Exception as e:
            _log(f"RDKit ligand: {e}")
    if not success:
        result.update({"status": "ok", "affinity": fpocket, "method": "Fpocket (ligand failed)"})
        return result

    # ── 2. 获取受体结构 ──────────────────────────
    rec = os.path.join(BASE, f"_dr_{_uid}.pdb")
    has_rec = False; source = ""; pdb_id = ""
    try:
        from uniprot_pdb import UNIPROT_PDB
        pdb_id = UNIPROT_PDB.get(uniprot, "")
    except ImportError: pass

    if receptor_pdb and len(receptor_pdb) > 3000:
        try:
            with open(rec, "w", encoding="utf-8") as f: f.write(receptor_pdb)
            has_rec = True; source = "frontend"
        except: pass
    if not has_rec and pdb_id and _download_pdb(pdb_id, rec):
        has_rec = True; source = f"RCSB/{pdb_id}"
    if not has_rec and _download_alphafold(uniprot, rec):
        has_rec = True; source = f"AlphaFold/{uniprot}"

    if not has_rec:
        try: os.remove(lig)
        except: pass
        result.update({"status": "ok", "affinity": fpocket, "method": "Fpocket (no structure)"})
        return result

    result["pdb"] = pdb_id or uniprot

    # ── 3. 确定对接中心：优先自动检测共晶配体 ──
    center = None
    # 自动检测：找含碳有机配体质心
    SKIP_RES = {"HOH","DMS","EDO","GOL","PEG","SO4","PO4","EPE","BME",
                "ACT","EOH","MPD","CIT","FMT","NO3","CL","NA","K","ZN",
                "MG","CA","MN","FE","CO","NI","CU"}
    all_het = {}
    try:
        with open(rec) as f:
            for line in f:
                if line.startswith("HETATM"):
                    resname = line[17:20].strip()
                    if resname in SKIP_RES: continue
                    try:
                        x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                        el = line[76:78].strip() or line[13:14].strip()
                        all_het.setdefault(resname, []).append((x,y,z,el))
                    except: pass
    except: pass
    best_c = None; best_n = 0
    for rn, atoms in all_het.items():
        if any(el=='C' for _,_,_,el in atoms) and len(atoms) > best_n:
            sx = sum(a[0] for a in atoms); sy = sum(a[1] for a in atoms)
            sz = sum(a[2] for a in atoms)
            best_c = (round(sx/len(atoms),1), round(sy/len(atoms),1), round(sz/len(atoms),1))
            best_n = len(atoms)
    # 计算蛋白中心（作为回退用）
    prot_cx = prot_cy = prot_cz = 0.0; prot_n = 0
    try:
        with open(rec) as f:
            for line in f:
                if line.startswith("ATOM"):
                    try: prot_cx+=float(line[30:38]); prot_cy+=float(line[38:46]); prot_cz+=float(line[46:54]); prot_n+=1
                    except: pass
    except: pass
    prot_center = (round(prot_cx/prot_n,1), round(prot_cy/prot_n,1), round(prot_cz/prot_n,1)) if prot_n>3 else (0,0,0)

    # ── 3. 确定口袋中心 ──
    # 优先级：手动标定 > 自动检测HETATM > 蛋白中心
    # force_blind 模式下跳过所有口袋检测，直接使用蛋白中心
    if force_blind:
        center = prot_center
        _log(f"center FORCE_BLIND=PROTEIN={center}")
    else:
        center = POCKET_CTR.get(pdb_id)
        if center:
            _log(f"center PRESET={center} pdb={pdb_id}")
        
        if not center and best_c:
            dx_m = best_c[0]-prot_center[0]; dy_m = best_c[1]-prot_center[1]; dz_m = best_c[2]-prot_center[2]
            dist_to_ctr = (dx_m**2+dy_m**2+dz_m**2)**0.5
            if dist_to_ctr > 15:
                _log(f"center HETATM={best_c} too far, using protein center")
                center = prot_center
            else:
                center = best_c
                _log(f"center HETATM={center}")
        
        if not center:
            center = prot_center
            _log(f"center PROTEIN={center}")

    # ── 4. 受体 PDBQT（完整受体，不截取）─────────
    recq = os.path.join(BASE, f"_drq_{_uid}.pdbqt")
    try:
        subprocess.run([OBABEL, rec, "-xr", "-opdbqt", "-O", recq],
                       timeout=120, capture_output=True)
        if not os.path.exists(recq) or os.path.getsize(recq) < 100:
            for f in [lig, rec, recq]: 
                try: os.remove(f)
                except: pass
            result.update({"status": "ok", "affinity": fpocket, "method": "Fpocket (bad rec)"})
            return result
    except subprocess.TimeoutExpired:
        for f in [lig, rec, recq]:
            try: os.remove(f)
            except: pass
        result.update({"status": "ok", "affinity": fpocket, "method": "Fpocket (rec timeout)"})
        return result

    # ── 5. Vina 对接 ─────────────────────────────
    out = os.path.join(BASE, f"_do_{_uid}.pdbqt")
    exhausted = 4 if source.startswith("AlphaFold") else 8
    # 大分子（>500Da）降低精度加速
    if len(smiles) > 80: exhausted = 2
    box = 30 if source.startswith("AlphaFold") else box_size
    
    def _run_vina(exh, sz):
        try:
            sp = subprocess.run([
                VINA, "--receptor", recq, "--ligand", lig, "--out", out,
                "--center_x", str(center[0]), "--center_y", str(center[1]), "--center_z", str(center[2]),
                "--size_x", str(sz), "--size_y", str(sz), "--size_z", str(sz),
                "--exhaustiveness", str(exh),
            ], capture_output=True, timeout=180, cwd=BASE)
            stdout = sp.stdout.decode('utf-8', errors='replace') if isinstance(sp.stdout, bytes) else (sp.stdout or '')
            stderr = sp.stderr.decode('utf-8', errors='replace') if isinstance(sp.stderr, bytes) else (sp.stderr or '')
            return stdout, stderr, None
        except subprocess.TimeoutExpired:
            return "", "", "timeout"

    _log(f"Vina try1 exhaust={exhausted} box={box} source={source}")
    stdout, stderr, err = _run_vina(exhausted, box)
    
    # 解析结合能
    affinity = None
    for line in stdout.split("\n"):
        ls = line.strip()
        if ls and ls[0].isdigit():
            parts = ls.split()
            try:
                if len(parts) >= 2 and 1 <= int(parts[0]) <= 20:
                    v = float(parts[1])
                    if -50 < v < 0: affinity = v; break
            except (ValueError, IndexError): continue
    
    # 第一次没找到，加大搜索框和算力重试
    if affinity is None and err != "timeout" and stdout:
        exh2 = min(exhausted * 2, 32)
        sz2 = min(box + 8, 30)
        _log(f"Vina retry exhaust={exh2} box={sz2}")
        stdout2, stderr2, err2 = _run_vina(exh2, sz2)
        for line in stdout2.split("\n"):
            ls = line.strip()
            if ls and ls[0].isdigit():
                parts = ls.split()
                try:
                    if len(parts) >= 2 and 1 <= int(parts[0]) <= 20:
                        v = float(parts[1])
                        if -50 < v < 0: affinity = v; break
                except (ValueError, IndexError): continue
        if err2 != "timeout":
            stdout = stdout2  # 用第二次的输出
    _log(f"parsed affinity={affinity}")

    # 若 RCSB Vina 失败，换 AlphaFold 再试
    if affinity is None and source.startswith("RCSB") and _download_alphafold(uniprot, os.path.join(BASE, f"_dr_af_{_uid}.pdb")):
        _log(f"Vina failed on RCSB, retry with AlphaFold for {uniprot}")
        af_recq = os.path.join(BASE, f"_drq_af_{_uid}.pdbqt")
        sp_prep = subprocess.run([OBABEL, os.path.join(BASE, f"_dr_af_{_uid}.pdb"), "-xr", "-opdbqt", "-O", af_recq], timeout=120, capture_output=True)
        if os.path.exists(af_recq) and os.path.getsize(af_recq) > 100:
            ax=ay=az=an=0
            with open(os.path.join(BASE, f"_dr_af_{_uid}.pdb")) as f:
                for line in f:
                    if line.startswith("ATOM"):
                        try: ax+=float(line[30:38]); ay+=float(line[38:46]); az+=float(line[46:54]); an+=1
                        except: pass
            c2 = (round(ax/an,1), round(ay/an,1), round(az/an,1)) if an>3 else center
            _log(f"AlphaFold center={c2}")
            sp_af = subprocess.run([
                VINA, "--receptor", af_recq, "--ligand", lig, "--out", out,
                "--center_x", str(c2[0]), "--center_y", str(c2[1]), "--center_z", str(c2[2]),
                "--size_x", "30", "--size_y", "30", "--size_z", "30",
                "--exhaustiveness", "8",
            ], capture_output=True, timeout=180, cwd=BASE)
            stdout_af = sp_af.stdout.decode('utf-8', errors='replace') if isinstance(sp_af.stdout, bytes) else (sp_af.stdout or '')
            if stdout_af:
                for line in stdout_af.split("\n"):
                    ls = line.strip()
                    if ls and ls[0].isdigit():
                        parts = ls.split()
                        try:
                            if len(parts) >= 2 and 1 <= int(parts[0]) <= 20:
                                v = float(parts[1])
                                if -50 < v < 0: affinity = v; break
                        except: continue
            if affinity is not None:
                source = f"AlphaFold/{uniprot}"
                center = c2
                stdout = stdout_af
                _log(f"AlphaFold retry OK: affinity={affinity}")

    # ── 7. 位姿分析 + 置信度评分 ────────────────
    pose_dist = None
    confidence = None
    if affinity is not None and os.path.exists(out) and os.path.getsize(out) > 100:
        try:
            pose_atoms = []
            with open(out) as f:
                in_model = False
                for line in f:
                    if line.startswith("MODEL"):
                        in_model = True; pose_atoms = []
                    elif line.startswith("ENDMDL"):
                        if in_model and pose_atoms: break
                        in_model = False
                    elif in_model and (line.startswith("ATOM") or line.startswith("HETATM")):
                        atom = line[12:16].strip()
                        if atom not in ("H", "HD", "HS", "HO"):
                            try:
                                pose_atoms.append((float(line[30:38]), float(line[38:46]), float(line[46:54])))
                            except: pass
            if pose_atoms:
                px = sum(a[0] for a in pose_atoms) / len(pose_atoms)
                py = sum(a[1] for a in pose_atoms) / len(pose_atoms)
                pz = sum(a[2] for a in pose_atoms) / len(pose_atoms)
                pose_dist = round(((px - center[0])**2 + (py - center[1])**2 + (pz - center[2])**2)**0.5, 1)
        except: pass
    if affinity is not None and pose_dist is not None:
        pose_penalty = max(0, 1.0 - pose_dist / 15.0)
        confidence = round(affinity * ps * pose_penalty, 3)
    _log(f"pose_dist={pose_dist} confidence={confidence} ps={ps}")

    # ── 8. 姿态 → SDF（仅有效负亲和力时） ──────
    pose_sdf = ""
    if affinity is not None and os.path.exists(out) and os.path.getsize(out) > 100:
        sdf_path = os.path.join(BASE, f"_do_{_uid}.sdf")
        try:
            subprocess.run([OBABEL, out, "-osdf", "-O", sdf_path],
                           timeout=30, capture_output=True)
            if os.path.exists(sdf_path) and os.path.getsize(sdf_path) > 100:
                with open(sdf_path) as f: pose_sdf = f.read()
            try: os.remove(sdf_path)
            except: pass
        except:
            try:
                with open(out) as f: pose_sdf = f.read()
            except: pass

    # ── 9. 受体可视化 ────────────────────────────
    # 直接发完整受体（不截取，避免 3Dmol 解析错误）
    receptor_frag = ""
    try:
        with open(rec) as f:
            receptor_frag = f.read()[:200000]
    except: pass

    method = f"Vina ({source})" if affinity is not None else f"Fpocket ({source})"
    _log(f"DONE method={method} pose_sdf_len={len(pose_sdf)}")

    # ── 10. 清理 ──────────────────────────────────
    for f in [lig, rec, recq, out]:
        try: os.remove(f)
        except: pass

    return {
        "status": "ok", "gene": gene,
        "affinity": affinity if affinity is not None else fpocket,
        "confidence": confidence,
        "pose_dist": pose_dist,
        "method": f"Vina ({source})" if affinity is not None else f"Fpocket ({source})",
        "pdb": pdb_id or uniprot,
        "pose_sdf": pose_sdf,
        "receptor_pdb": receptor_frag,
        "pocket_center": list(center),
        "ligand_smiles": smiles,
        "fpocket": fpocket,
        "pocket_score": ps,
    }
