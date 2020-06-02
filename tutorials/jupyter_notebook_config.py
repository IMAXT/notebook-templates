# Reference: https://svds.com/jupyter-notebook-best-practices-for-data-science/
from subprocess import check_call
from pathlib import Path

GIT = False

def post_save(model, os_path, contents_manager):
    """post-save hook for converting notebooks to .py scripts"""
    if model['type'] != 'notebook':
        return # only do this for notebooks
    if os_path.endswith('.py'):
        return
    if 'Untitled' in os_path:
        return
    p = Path(os_path)
    for j in range(9):
        try:
            f = p.parent / f'{p.stem}.py~{3-j}~'
            f.rename(p.parent / f'{p.stem}.py~{4-j}~')
        except:
            pass
    try:
        f = p.parent / f'{p.stem}.py'
        f.rename(p.parent / f'{p.stem}.py~{1}~')
    except:
        pass
    try:
        check_call(['jupyter', 'nbconvert', '--to', 'script', p.name], cwd=str(p.parent))
        if GIT:
            try:
                check_call(['git', 'add', f'{p.stem}.py'], cwd=str(p.parent))
                check_call(['git', 'commit', '-m', f'"Update {p.name}"'])
            except:
                pass
    except:
        pass

c.FileContentsManager.post_save_hook = post_save

