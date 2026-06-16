# -*- coding: utf-8 -*-
"""
TCM-TargetMiner 数据库修复工具
删除旧数据库，重新初始化
"""

import os

BASE = os.path.dirname(os.path.realpath(__file__))
DB_PATH = os.path.join(BASE, "tcm_demo.db")

print("=" * 50)
print("  TCM-TargetMiner 数据库修复工具")
print("=" * 50)

# 1. 删除旧数据库
if os.path.exists(DB_PATH):
    size = os.path.getsize(DB_PATH)
    print(f"  发现旧数据库: {size} 字节")
    confirm = input("  是否删除重建？(y/n): ")
    if confirm.lower() == 'y':
        os.remove(DB_PATH)
        print("  ✓ 旧数据库已删除")
    else:
        print("  ✗ 操作取消")
        exit()
else:
    print("  未找到数据库文件，将新建")

# 2. 运行 app.py 的 init_db 重建数据库
print("\n  正在导入数据...")
# 直接导入 app 模块会启动 Flask，所以我们只执行 init_db
import sqlite3
import json
import math
import sys

# 加载 app 模块但不启动 Flask
os.chdir(BASE)
sys.path.insert(0, BASE)

# 手动执行 init_db 的逻辑
from app import init_db, get_db

init_db()

# 验证
conn = get_db()
c = conn.cursor()
c.execute("SELECT COUNT(*) FROM herbs")
herb_count = c.fetchone()[0]
c.execute("SELECT COUNT(*) FROM ingredients")
ing_count = c.fetchone()[0]
c.execute("SELECT COUNT(*) FROM targets")
tgt_count = c.fetchone()[0]
conn.close()

print(f"\n  ✓ 数据库重建完成！")
print(f"    中药: {herb_count} 味")
print(f"    成分: {ing_count} 个")
print(f"    靶点: {tgt_count} 个")
print(f"\n  现在可以启动服务器了：")
print(f"    双击「一键启动.bat」")
print(f"    或运行: python app.py")
print("=" * 50)
input("\n按回车键退出...")
