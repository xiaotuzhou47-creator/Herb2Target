import urllib.request
import os

files = [
    ('https://unpkg.com/vue@3/dist/vue.global.prod.js', 'vue3.js'),
    ('https://unpkg.com/element-plus/dist/index.css', 'element-plus.css'),
    ('https://unpkg.com/element-plus/dist/index.full.min.js', 'element-plus.js'),
    ('https://unpkg.com/@element-plus/icons-vue/dist/index.iife.min.js', 'element-icons.js'),
    ('https://unpkg.com/axios/dist/axios.min.js', 'axios.js'),
]

dst = os.path.dirname(os.path.abspath(__file__))
os.makedirs(dst, exist_ok=True)

for url, fname in files:
    local = os.path.join(dst, fname)
    try:
        urllib.request.urlretrieve(url, local)
        print(f'OK: {fname}')
    except Exception as e:
        print(f'FAIL: {fname} -> {e}')
