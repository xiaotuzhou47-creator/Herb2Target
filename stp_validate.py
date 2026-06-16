# -*- coding: utf-8 -*-
"""SwissTargetPrediction validation via browser automation"""
import sys, json, time, re, os
sys.stdout.reconfigure(encoding='utf-8')

from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.options import Options
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC

VALIDATION = [
    ("Resveratrol", "C1=CC(=CC=C1C=CC2=CC(=CC(=C2)O)O)O", "PTGS2"),
    ("Curcumin", "COC1=C(C=CC(=C1)C=CC(=O)CC(=O)C=CC2=CC(=C(C=C2)O)OC)O", "TNF"),
    ("Berberine", "COC1=C(C2=C(C=C1)C=C3C4=CC5=C(C=C4CC[N+]3=C2)OCO5)OC", "AKT1"),
    ("Baicalein", "C1=CC=C(C=C1)C2=CC(=O)C3=C(C=C(C(=C3O2)O)O)O", "PTGS2"),
    ("Kaempferol", "C1=CC(=CC=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)O)O", "AKT1"),
    ("Hesperetin", "COC1=C(C=CC(=C1)C2CC(=O)C3=C(C=C(C=C3O2)O)O)O", "PTGS2"),
    ("Apigenin", "C1=CC(=CC=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O", "PTGS2"),
    ("Wogonin", "COC1=C(C=C(C2=C1OC(=CC2=O)C3=CC=CC=C3)O)O", "PTGS2"),
    ("EGCG", "C1=C(C=C(C(=C1O)O)O)C2CC(=O)C3=C(C=C(C(=C3O2)O)O)O", "EGFR"),
    ("Luteolin", "C1=CC(=C(C=C1C2=CC(=O)C3=C(C=C(C=C3O2)O)O)O)O", "TNF"),
]

OUTPUT_FILE = "E:/Herb2Target_manuscript/stp_validation_20260615.json"

options = Options()
options.binary_location = r"C:\Program Files\Google\Chrome\Application\chrome.exe"
options.add_argument("--no-sandbox")
options.add_argument("--disable-dev-shm-usage")
options.add_argument("--window-size=1200,900")

driver = webdriver.Chrome(options=options)
wait = WebDriverWait(driver, 60)

results = []

try:
    for idx, (name, smiles, target) in enumerate(VALIDATION):
        print(f"\n[{idx+1}/10] {name} -> {target}")
        print(f"  SMILES: {smiles}")

        driver.get("https://www.swisstargetprediction.ch/")
        time.sleep(2)

        try:
            box = wait.until(EC.presence_of_element_located((By.ID, "smilesBox")))
            box.clear()
            box.send_keys(smiles)
            time.sleep(1)
        except Exception as e:
            print(f"  ERROR: Cannot find SMILES input: {e}")
            results.append({"compound": name, "expected_target": target, "stp_rank": None, "error": str(e)})
            continue

        try:
            form = driver.find_element(By.ID, "myForm")
            driver.execute_script("document.getElementById('myForm').submit();")
        except Exception as e:
            print(f"  ERROR: Cannot submit: {e}")
            results.append({"compound": name, "expected_target": target, "stp_rank": None, "error": str(e)})
            continue

        print("  Waiting for prediction...")
        time.sleep(15)

        driver.save_screenshot(f"E:/demo/demo/stp_{name.lower()}.png")
        html = driver.page_source
        with open(f"E:/demo/demo/stp_{name.lower()}.html", "w", encoding="utf-8") as f:
            f.write(html)

        # Try to find rank of expected target from page text
        page_text = driver.find_element(By.TAG_NAME, "body").text

        # Look for expected target in the text
        target_upper = target.upper()
        lines = page_text.split("\n")
        found_idx = None
        for li, line in enumerate(lines):
            if target_upper in line.upper() or target_upper in line:
                found_idx = li + 1
                print(f"  Found '{target}' in line {li}: {line[:100]}")
                break

        # Also try looking for table structure
        rows = driver.find_elements(By.XPATH, "//tr")
        for ri, row in enumerate(rows):
            cells = row.find_elements(By.TAG_NAME, "td")
            cell_text = " ".join(c.text.strip().upper() for c in cells)
            if target_upper in cell_text:
                found_idx = ri + 1
                break

        results.append({"compound": name, "expected_target": target, "stp_rank": found_idx})

        if found_idx:
            print(f"  >>> Rank: #{found_idx}")
        else:
            print(f"  >>> Target not found in results (may be >15)")
            results[-1]["stp_rank"] = ">15"

        time.sleep(2)

finally:
    driver.quit()

# Summary
print("\n" + "=" * 50)
print("STP VALIDATION RESULTS - 2026-06-15")
print("=" * 50)

rank_values = []
for r in results:
    if r["stp_rank"] == ">15":
        rank_values.append(16)
    elif r["stp_rank"] is not None:
        rank_values.append(r["stp_rank"])
    else:
        rank_values.append(16)

top5 = sum(1 for v in rank_values if v <= 5)
top10 = sum(1 for v in rank_values if v <= 10)

for r in results:
    rank_str = str(r["stp_rank"]) if r["stp_rank"] else "N/A"
    print(f"  {r['compound']:15s} -> {r['expected_target']:6s}: Rank {rank_str}")

print(f"\nSTP Top-5:  {top5}/10 ({top5*10}%)")
print(f"STP Top-10: {top10}/10 ({top10*10}%)")

output = {
    "date": "2026-06-15",
    "tool": "SwissTargetPrediction",
    "results": results,
    "summary": {"top5": f"{top5}/10", "top10": f"{top10}/10"}
}
with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
    json.dump(output, f, ensure_ascii=False, indent=2)
print(f"\nSaved to: {OUTPUT_FILE}")
