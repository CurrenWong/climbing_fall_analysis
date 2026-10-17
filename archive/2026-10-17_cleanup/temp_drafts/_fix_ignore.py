import io
import os
import subprocess

BASE = r'D:\Project\climbing_fall_analysis'
gi = os.path.join(BASE, '.gitignore')
txt = io.open(gi, encoding='utf-8').read() if os.path.exists(gi) else ''
if '*_TODO.md' not in txt:
    io.open(gi, 'a', encoding='utf-8').write(
        '\n# 会话级草稿（Coding Mode 约定：不入库）\n*_TODO.md\n')
    print('gitignore += *_TODO.md')
else:
    print('gitignore 已含 *_TODO.md')

files = subprocess.run(['git', 'ls-files', '*_TODO.md'], cwd=BASE,
                       capture_output=True, text=True).stdout.split()
root = [f for f in files if '/' not in f and '\\' not in f]
print('已跟踪 TODO 总数:', len(files), '｜根目录:', root)
if root:
    subprocess.run(['git', 'rm', '--cached', '-q'] + root, cwd=BASE)
    print('已 git rm --cached', len(root), '个（仅根目录）')
